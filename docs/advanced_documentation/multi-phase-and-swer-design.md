<!--
SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>

SPDX-License-Identifier: MPL-2.0
-->

# Design proposal: reduced-phase components and SWER networks

```{note}
This is a design and investigation document. It describes what is needed to support one- and two-phase
components in power flow and short-circuit calculations. Milestones M0 to M11 are implemented. One- and two-phase
components and delta loads stay an experimental feature (decisions D9 and D45). For how to use the feature, see
[Reduced-phase networks](../user_manual/reduced-phase-networks.md).
```

## Context

Many distribution networks are not fully three-phase. Rural networks in Australia are a typical example:

* Two-phase and single-phase laterals branch off three-phase feeders.
* **SWER** (single-wire earth return) systems are fed from an isolating transformer, phase-to-phase on the
  primary (e.g. 22 kV a-b) and phase-to-earth on the secondary (12.7 kV or 19.1 kV).
* The SWER line is one conductor; the earth is the return path.
* Customers are supplied by single-phase distribution transformers with a 240 V or 240-0-240 V (split-phase)
  secondary.
* Loads can be phase-to-earth or phase-to-phase.

Power-grid-model (PGM) can model none of this today. This document:

1. lists where the code base assumes three phases;
2. proposes a numerical strategy and data model to remove those assumptions for power flow and short circuit;
3. breaks the work into milestones.

Design decisions taken as input:

* **Scope:** the full SWER chain (lines, isolating and distribution transformers, split-phase secondaries and
  phase-to-phase loads).
* **Data model:** extend the existing `asym_line` with a `phases` attribute instead of adding a new line type.
  Node phases are derived from what is connected to them, so nodes get no new attribute.

In this document `<core>` means `power_grid_model_c/power_grid_model/include/power_grid_model`.

## Current architecture and its three-phase assumptions

### Three phases are fixed in the types

* **Core types.** `ComplexValue<asymmetric_t>` is an Eigen `Array<T, 3, 1>` and `ComplexTensor<asymmetric_t>` is
  an `Array<T, 3, 3>` (`<core>/common/three_phase_tensor.hpp:30-33,129-143`). `static_assert`s fix their memory
  layout (lines 146-161).
* **Block sizes.** The sparse block size is `is_symmetric_v<sym> ? 1 : 3` (`<core>/math_solver/block_matrix.hpp:20`).
  The asymmetric Y-bus therefore uses 3×3 complex blocks, and the Newton-Raphson Jacobian (`PFJacBlock`,
  `<core>/math_solver/newton_raphson_pf_solver.hpp:206-225`) uses 6×6 real blocks.
* **Data interface.** Asymmetric attributes have type `CType::c_double3` (`<core>/common/enum.hpp:121`). The C API
  uses `PGM_double3` and Python uses `(3,)f8` (`src/power_grid_model/_core/power_grid_meta.py:42-45`).
* **Phase selection.** No component has a per-phase connection mask. The only phase selector in the model is
  `fault.fault_phase`, and it only affects which rows the short-circuit solver modifies.
* **Status.** Branch `from_status`/`to_status` and appliance `status` switch all three phases together.

### Missing phases break the numerics

* **Singular blocks.** A phase with no admittance gives a zero row and column in a Y-bus diagonal block. The block
  LU factorisation then throws `SparseMatrixError` (`<core>/math_solver/sparse_lu_solver.hpp:160-163`). Power flow
  and short circuit never enable pivot perturbation; only state estimation does.
* **Zero voltages.** A phase voltage of zero causes several problems:
  * `conj(s / u)` returns NaN in the iterative current solver
    (`<core>/math_solver/iterative_current_pf_solver.hpp:180-190`) and in result calculation
    (`<core>/math_solver/common_solver_functions.hpp:153,407`).
  * `calculate_hnml` scales by $V_i V_j$ (`newton_raphson_pf_solver.hpp:498-507`), so the Newton-Raphson Jacobian
    rows for that phase become zero.
  * The multiplicative voltage update keeps $V = 0$ for ever.

### Component models assume balanced or full three-phase construction

* **`asym_line`** (`<core>/component/asym_line.hpp`) is the only branch built natively in abc.
  * It always inverts the full 3×3 series impedance (`y_series_abc_ = inv(z_series)`, line 39). A missing phase
    makes this inverse singular or NaN.
  * The open-end shunt is only applied when *every* element of $Y_{shunt}$ is non-zero (line 110 at the time;
    now `asym_line.hpp:269`). Known pre-existing behaviour; Michael decided (2026-10-15) to leave it as it is for
    now.
  * The neutral is removed by Kron reduction (`<core>/component/line_utils.hpp:12-19`), which assumes the neutral
    is at zero volts.
* **Sequence-domain components** (line, link, shunt, source and transformer) build their abc matrices from
  positive- and zero-sequence values:
  * line, link and shunt: `<core>/component/edge.hpp:227-245`;
  * source: `<core>/calculation_parameters.hpp:220-231`;
  * transformer: `<core>/component/transformer.hpp:333-340`, which also enforces three-phase winding types.
* **Initial voltages** come from the rotating constructor `Vector(T x)` → `(x, x·a², x·a)`
  (`three_phase_tensor.hpp:47-53`), plus one scalar `phase_shift` per bus (`MathModelTopology`,
  `<core>/calculation_parameters.hpp:165-187`). That rotation is wrong for a split-phase bus, whose two legs are
  180° apart.
* **Loads** are wye (phase-to-earth) only (`<core>/component/load_gen.hpp:124-139`). There is no delta connection.
* **Sources** are balanced three-phase Norton equivalents with a scalar `u_ref`.

### Short circuit

* **Symmetry choice.** The symmetric solver is used only when every fault is `three_phase`
  (`<core>/main_model_impl.hpp:362-368`).
* **Phase indices.** Fault phases map to fixed indices 0–2 (`set_phase_index`,
  `<core>/math_solver/short_circuit_solver.hpp:392-433`). Bolted faults substitute unknowns (lines 144-207), and
  faults with impedance are added to the diagonal block (lines 209-253).
* **One fault type per calculation.** All faults in one calculation must share the same type and phase
  (`check_input_valid`, lines 444-467).

### Pre-existing defects found during this investigation

* **Symmetric parameters of `asym_line`.** `average_of_off_diagonal_of_matrix`
  (`<core>/common/matrix_utils.hpp:16-18`) adds `matrix(1, 2)` twice and leaves out `matrix(0, 2)`.
  `AsymLine::sym_calc_param` uses it, so symmetric results for any non-transposed `asym_line` are slightly wrong.
  This is fixed in milestone M0.
* **Sign of the mutual capacitance in the Python validation (decision D40).** The rule
  `c_ba`, `c_ca`, `c_cb` $\geq 0$ contradicts the model of the core, which uses the c matrix as the nodal (Maxwell)
  matrix: `c0` and `c1` give an off-diagonal of $(c_0 - c_1)/3 < 0$, and every OpenDSS-validated case has negative
  mutual capacitances. M9 drops the sign rule, in a commit of its own so that it can be proposed upstream separately.

### What is generated and what is written by hand

* **Generated.** Attribute structs, metadata, the C API dataset definitions and the Python
  `dataset_definitions.py` all come from `code_generation/data/attribute_classes/*.json` and
  `code_generation/data/dataset_class_maps/dataset_definitions.json` via `code_generation/code_gen.py`.
* **Written by hand:**
  * enums (`<core>/common/enum.hpp`, `src/power_grid_model/_core/enum.py`);
  * component classes, plus `<core>/all_components.hpp` and `<core>/main_core/main_model_type.hpp`;
  * the math layer;
  * Python validation (`src/power_grid_model/validation/_validation.py`);
  * docs and tests.
* **Template.** The PR that introduced `asym_line` (merge commit `9cf8afb98`) is a good model for adding a
  component.

## Recommended numerical strategy: phase mask with decoupled identity equations

**Keep the fixed 3×3 asymmetric block structure.** Each bus keeps three phase slots (a, b, c). A slot with no
conductor is *absent*. Each absent slot gets a decoupled trivial equation $1 \cdot U_p = 0$. The solvers and the
sparse LU never see a singular block, and ordinary three-phase grids are unaffected.

### Phase mask

* **Phase mask type.** `PhaseMask` is a one-byte wrapper class in `<core>/calculation_parameters.hpp`: bit 0 for a, bit 1
  for b, bit 2 for c. All bit operations live inside the class.
* **New fields.** Added to `MathModelTopology`:
  * `std::vector<PhaseMask> bus_phase_mask`: one entry per bus.
  * `bool has_reduced_phase`: a fast-path flag. When it is false the mask stays empty and nothing below runs.
    Read the mask through `bus_phase(bus)`, which returns abc when the flag is false.
  * `std::vector<SlotPattern> bus_slot_pattern` (see [Initial voltages](#initial-voltages)) is not part of the
    topology milestone; it was added together with the start values in M6 (decision D25). It is filled exactly when
    `bus_phase_mask` is, and read through `slot_pattern(bus)`, which returns `three_phase` when the flag is false.
* **Incidence per branch.** A virtual `Edge::phase_incidence()` returns a `BranchPhaseIncidence`: the non-earth
  terminals at the from and to side, plus the coupling. `per_phase` couples slot k on one side to slot k on the other
  (line-like, both sides equal). `winding` couples all terminals of a side through a winding (transformer-like). It
  defaults to abc → abc, per phase; `Transformer` and `GenericBranch` report a three-phase winding, and each leg of a
  three-winding transformer is a three-phase winding. The values are registered in
  `ComponentConnections::branch_phase_incidence` (`<core>/main_core/topology.hpp`), which stays empty when every
  branch is abc → abc.
* **Propagation.** `Topology::propagate_phases()` runs in `Topology::build_topology()` after `dfs_search()` and returns
  immediately when `branch_phase_incidence` is empty. Otherwise:
  1. Build a directed phase graph with the same edges as the global graph.
  2. Seed *every* connected source bus (not only the slack bus) with abc, then propagate to a fixed point with a
     worklist. Masks only grow, so this terminates.
     * **Line-type branches:** a present slot on one side makes the same slot present on the other side.
     * **Transformer windings:** the far side becomes present only if *all* non-earth terminals of the feeding
       winding are present. Otherwise a winding could back-feed an absent slot and create a floating, singular
       sub-problem.
     * **Links and supernodes:** on the v2 path, nodes are merged before the topology is built, so a merged node
       carries the union. On the legacy path, a link reports the coupling `PhaseCoupling::link` (abc → abc). It
       propagates like a per-phase branch, which gives the same union, but it is an ideal connection that carries
       whatever phases its nodes have, so it is never validated as dangling (decision D2).
  3. Validate after the fixed point (`check_branch_phases`), so that the outcome does not depend on the visit order:
     * **Windings.** A connected winding side (also of a branch that is connected on one side only) must have *all*
       of its terminals present at its node. Otherwise it raises `MissingPhaseForWinding`. A winding with none of its
       terminals present is an error as well (decision D1): it would stamp its admittance on absent slots.
     * **Dangling conductors.** At every connected side of a per-phase branch whose node is in a math model, the
       branch terminals must be a subset of the node phases. Otherwise it raises `DanglingConductor` with the branch
       index, the dangling phases and the phases at both ends. This covers a three-phase line between two separately
       fed single-phase nodes, a `b` line off an `a` node, and an abc stub hanging from an `a` node. For a per-phase
       branch connected on both sides, propagation already makes the present terminals equal at both ends, so a
       conductor is either present at both ends or dangling at both ends.
     * **Links** are skipped.
  4. Write the masks. A math model gets `has_reduced_phase` and a mask vector only if at least one of its buses is not
     abc. With the two rules above, every bus in a math model is reached by at least one phase from a source, so a bus
     with the mask `none` cannot occur.
* **Rebuilds.** `phases` and the connection attributes are input-only. The mask therefore changes only through
  status updates, which already trigger a topology rebuild.

### Solver changes

| Area | Change |
|------|--------|
| Y-bus (`<core>/math_solver/y_bus.hpp`, `update_admittance_entries`) | Asymmetric only. **Mask, then identity:** every entry $(i, j)$ becomes $M_i Y_{ij} M_j$, and on the diagonal $(p, p)$ is *overwritten* with 1 for each absent `p`. Adding 1.0 would only be correct if every contribution were exactly zero on absent slots, which sequence-based shunts (mutual terms) and legacy links (1e6 on the diagonal) break. With masking, each absent row is exactly $e_p$ whatever the component data. The incremental update path reuses the same function. Branch and shunt flows are masked by the phases of the bus they flow into. Masking is by assignment, never by multiplication, so a NaN in an absent slot is discarded too. The identity rows are part of `YBus::admittance()` (decision D8). |
| Newton-Raphson (`newton_raphson_pf_solver.hpp`) | New `apply_absent_phase_constraints()`, called after `apply_pv_constraints()`. For each absent slot it zeroes the P and Q rows in all entries of the bus row, sets $H(p,p) = L(p,p) = 1$ and sets the mismatch to 0. The update keeps $V_p = 0$ and never divides by it. A PV bus (voltage regulator) on a bus without all phases raises `ReducedPhaseNotSupported`. |
| Newton-Raphson starts | The flat, average-source and linear starts mask the initial voltage of absent slots to exactly 0. The linear start also masks the load admittance per bus. Since M6 the flat and average-source starts first apply the slot patterns (`bus_slot_pattern`, decision D25), so that slot b of a split-phase bus starts opposite to slot a. |
| Iterative current (`iterative_current_pf_solver.hpp`) | Apply the slot patterns (decision D25) to the average-source start, then mask the initial voltage and, after the load and source loop, the right-hand side (which holds `conj(s/0)` = NaN on absent slots). The prefactorised matrix already has the identity rows. `linear_current` is covered as well. |
| Linear power flow (`common_solver_functions.hpp`, `add_linear_loads`) | The load admittance $-\overline{s}$ is masked per bus. Otherwise it would be added to the identity on absent slots and give $(1 - \overline{s_p}) U_p = 0$, which is singular for $s_p = 1$. |
| Newton-Raphson, delta loads (M7) | `add_loads` dispatches a delta load_gen (`topo.is_delta_load_gen`) to `add_delta_load`, which adds the branch powers and the full Jacobian (angle and magnitude terms, phases coupled within the bus's diagonal block). The linear start adds the delta constant admittance. See "Phase-to-phase loads" below. |
| Iterative current, delta loads (M7) | The right-hand side gets $D^T I_\Delta(U)$ instead of the wye current; the matrix is unchanged. |
| Linear power flow, delta loads (M7) | `add_linear_loads` adds the delta constant admittance $-D^T \mathrm{diag}(\overline{S_k}/U_{\Delta,nom}^2) D$ to the bus's diagonal block. |
| Results (`common_solver_functions.hpp`) | The load/gen `s` and `i` (`conj(s/u)`) are set to 0 on absent slots. Bus injection is exactly 0 there because $U_p = 0$. For a delta load_gen (M7), `calculate_load_gen_result` gives the phase current $D^T I_\Delta$ and $s = U \circ \overline{I}$ (decision D30). |
| Short circuit (`short_circuit_solver.hpp`) | The Y-bus identity carries through `copy_y_bus`. Fault substitution only touches fault phases. In a reduced-phase model the solver checks that the fault phases are present at the faulted bus and raises `InvalidShortCircuitPhases` otherwise, so a three-phase fault, or a fault on an absent phase, at a reduced bus is rejected instead of silently giving zero fault current. Since M8 the model-level check with the fault ID (D34) runs first; this is the backstop for direct solver use. |
| State estimation (`math_solver.hpp`) | Not supported at first: `MathSolver::run_state_estimation` raises `ReducedPhaseNotSupported` when `has_reduced_phase`, for both symmetries and all methods. The observability check works per bus, not per phase, and pivot perturbation would hide modelling errors. |
| Symmetric calculations (`math_solver.hpp`) | The positive-sequence Y-bus cannot represent absent phases. `MathSolver<symmetric_t>` raises `ReducedPhaseNotSupported` for power flow and short circuit on a reduced-phase model. This is a permanent backstop behind the model-level guard of M4. |
| Voltage regulators / PV buses | Not supported on reduced-phase buses at first (see Newton-Raphson). |
| Outputs (`<core>/main_core/output.hpp`, `<core>/component/node.hpp`) | No change. The solvers leave every absent slot at exactly 0, so the node reports 0 for `u`, `u_pu`, `u_angle` ($\arg 0 = 0$), `p` and `q`, and every flow (branch, shunt, load/gen) is 0. See [Absent-phase output](#absent-phase-output). |

### Alternatives considered

* **A tiny regularising admittance to earth on absent slots.** Rejected:
  * The voltages it produces are not physically meaningful.
  * The Newton-Raphson Jacobian becomes ill-conditioned.
  * Outputs need a tolerance-based test to decide what is absent.
  * Loads on absent phases would silently draw current.
* **Variable block size or per-phase scalar unknowns.** Rejected. Every `ComplexValue<asymmetric_t>` user, the
  `Block` and `SparseLUSolver` templates, `YBusStructure` and `PFJacBlock` would need a run-time size. That is
  effectively a rewrite of the asymmetric math, and normal three-phase grids would get slower. The saving is only
  about one third of the unknowns on SWER buses.

The identity approach costs one trivial pivot per absent slot. It costs nothing when `has_reduced_phase` is false.

### Absent-phase output

**Policy: 0 everywhere.** An absent phase slot reports exactly 0 for every output field: node `u`, `u_pu`,
`u_angle`, `p` and `q`, and all branch, shunt and appliance flows. This is the same convention as for a
de-energised component (`Node::get_null_output`). In PGM, NaN means "input not provided", not "no conductor".

* **Why not NaN for absent voltages.** NaN would keep "absent" apart from "present at 0 V", which matters most in
  short-circuit results. But it would make NaN mean something new in output, need NaN-aware comparison in the
  Python and C++ validation tooling (per-element `equal_nan`, a non-NaN angle reference, a per-phase instead of
  per-element skip), and force users who aggregate voltages to use NaN-aware min/max.
* **How to tell absent from 0 V.** The node phases are derived from the input (`asym_line.phases` and the
  transformer windings), not reported. If users need the distinction, a future output attribute could expose the
  derived node phases.
* **Implementation.** No output code is needed. The solvers keep absent slots at exactly 0 (a positive zero, so
  $\arg 0 = 0$ and never $\pi$), and the node and branch output formulas map 0 to 0. Flows stay physically exact,
  so `loading`, `pf` and per-phase sums stay finite.

### Decisions

#### Taken in M2

Absent phases and the first guards:

* **D1. Winding without terminals at a connected side is an error.** This overturns the M1 rule that allowed it.
  Otherwise the winding would stamp its admittance on absent slots, and a bus with mask `none`, whose loads are
  dropped silently, could appear.
* **D2. `PhaseCoupling::link`.** Legacy-path links are phase-agnostic and exempt from the dangling check. Without
  it, a link between two SWER nodes on the legacy path (models with node-injection sensors) would raise
  `DanglingConductor`.
* **D5. Absent-phase output** is 0 everywhere (see above), not the NaN-for-voltages hybrid that was proposed.
* **D6. One exception class.** `ReducedPhaseNotSupported` is used for state estimation, symmetric calculations and
  voltage control on reduced-phase buses. M4 reuses it instead of adding `SymmetricCalculationNotSupported`. M4 also
  adds the Python error mapping (`power_grid_model.errors`) for `ReducedPhaseNotSupported`, `DanglingConductor` and
  `MissingPhaseForWinding`.

#### Taken in M4

Routing:

* **D7. Order of M3 and M4.** M4 lands before M3, so that reduced-phase models are never routed to the symmetric
  solver. A model with reduced-phase components always runs the asymmetric short circuit, also when every fault is
  `three_phase` and when there is no fault at all.

#### Taken in M3

`asym_line` with `phases`, the Y-bus and the experimental gate. Decided by recommendation; the points marked "decided
by Michael" were reviewed on 2026-10-15:

* **D8. Identity rows in the Y-bus.** `YBus::admittance()` keeps the identity rows of absent slots; this is
  documented at the accessor. Consumers that need the physical admittance (future state estimation) must mask these
  slots.
* **D9. Experimental gate** (decided by Michael). Until the maintainers graduate the feature (decision D45), a
  calculation on a model with one- or two-phase branches raises `ExperimentalFeature` ("A network with
  one- or two-phase branches is an experimental feature") unless experimental features are enabled. The gate is
  checked by the C API before the calculation starts, so it runs *before* the M4 guard: with experimental features
  disabled, every calculation type gives `ExperimentalFeature`; with them enabled, the M4 rules apply
  (`ReducedPhaseNotSupported` for symmetric power flow and state estimation). Model construction, including the
  `phases` validation below, is not gated, so invalid input is reported regardless.
* **D10. `PhaseConnection` numbering** (decided by Michael). The values are aligned
  with `FaultPhase`: `abc = 0, a = 1, b = 2, c = 3, ab = 4, ac = 5, bc = 6`. Not provided (`na`) means `abc`; there
  is no `nan` enumerator.
* **D11. Validation of `asym_line` input with `phases`** (decided by Michael).
  An entry of the r, x or c matrix is required if both of its conductors are present and must be empty otherwise.
  A non-empty value on an absent phase is rejected rather than ignored, so that `phases = a` set by mistake on a
  three-phase line does not silently drop two phases. The present c entries are required (0 is allowed), as for
  `abc`. `c0`/`c1` are only allowed with `abc`. Errors raise `InvalidPhaseConnection` with the line ID. An unknown
  `phases` value raises `MissingCaseForEnumError`.
* **IDs in topology errors.** `DanglingConductor` and `MissingPhaseForWinding` name the branch by ID ("branch with
  ID 42", "branch3 with ID 7") when the model registered the IDs, which it does whenever a branch is not abc → abc.
  Hand-built topologies without IDs keep the index ("branch #2").

#### Taken in M3, for loads and shunts

Appliances at nodes without all phases:

* **D3. Load/gen power on an absent phase is a model-level error** (`LoadGenOnAbsentPhase`, with the component ID).
  * **Where:** in the asymmetric power-flow input preparation
    (`<core>/main_core/calculation_input_preparation.hpp`), only when a math model has `has_reduced_phase`.
  * **What counts:** per absent phase, the real or imaginary part of the injection the calculation uses is not NaN
    and not 0 (an exact test). NaN and 0 both mean "no power". A component with status off has no injection, and a
    de-energised component is skipped.
  * **`sym_load`/`sym_gen` on a node without all phases is rejected** when its power is not 0 (decided by
    Michael). A `sym_load` is a balanced three-phase load: in an asymmetric
    calculation it puts $s/3$ on every phase. The alternatives (ignore the absent phases, or spread the power over
    the present phases) would silently drop load, or make one `sym_load` change meaning between scenarios when the
    node phases change with a status update. Rejecting first is forward-compatible. The message tells the user to
    use an `asym_load`/`asym_gen` with zero power on the absent phases.
  * **Batch:** the check runs per scenario on the updated values. The batch cache run swallows
    `LoadGenOnAbsentPhase` (decided by Michael), like `SparseMatrixError`, because
    the update data may fix the value; a value that no scenario fixes then fails each scenario on its own.
* **D4. A shunt with status on at a node without all phases raises `ShuntAtReducedPhaseNode`** ("Shunt {id} is
  connected to a node without all phases. A shunt at a node without all phases is not supported."), in the asymmetric
  power-flow and the short-circuit input preparation, until an asymmetric shunt exists. A shunt with status off is
  allowed. The Y-bus masking of M2 stays as a defensive layer. The error has its own class, so that the batch cache
  run can swallow it like `LoadGenOnAbsentPhase` without swallowing `ReducedPhaseNotSupported`. Each scenario is then
  judged on its own: a base model with an energised shunt at a reduced-phase node succeeds in a batch where every
  scenario switches the shunt off (decided by Michael: consistent with loads; revisit if a use case needs
  otherwise).

#### Taken in M5

The `single_phase_transformer` (all decided by Michael, 2026-10-15):

* **D12. `WindingConnection1P` values.** `a = 1, b = 2, c = 3, ab = 4, ac = 5, bc = 6`, the numbers of
  `PhaseConnection`/`FaultPhase` without `abc`. The connection is required: not provided (`na`), 0 and unknown values
  raise `MissingCaseForEnumError`. The first letter is the positive terminal: `ab` is $U_a - U_b$, `ac` is
  $U_a - U_c$ (the design's `ca` is `ac` with reversed polarity), and `a` is $U_a - U_E$ with $E$ the earthed end.
  The centre-tapped values are not part of M5; M6 appends only `ab_center_tapped = 7` (decision D19), which is
  ABI-compatible.
* **D13. Polarity through `clock`.** A `clock` attribute (`IntS`, not provided means 0) expresses the relative
  polarity: any value equal to 0 or 6 modulo 12 is accepted, 6 reverses the to winding; any other value raises
  `InvalidTransformerClock`. This covers phase-to-earth and phase-to-phase windings uniformly and keeps the enum
  set-like.
* **D14. Earthing.** On a phase-to-earth side the earthed terminal is connected to earth through
  `r/x_grounding_<side>`; not provided means 0, solidly earthed (as for `transformer`). A phase-to-phase side has no
  earthed terminal: a non-empty `r/x_grounding_<side>` raises `InvalidPhaseConnection` with the ID (following D11,
  reject rather than ignore). Every connected phase-to-phase side gets the low-susceptance fallback: a common-mode
  admittance $\tfrac{y_{lo}}{2}\mathbf{1}_{2\times2}$ on its two slots, in per unit, with
  $y_{lo} = -j \cdot 10^{-8} \cdot s_n / S_{b,1p} / u_k$ (the two-terminal analogue of the zero-sequence
  $y_0/3 \cdot \mathbf{1}_{3\times3}$ of `transformer`); the differential mode is untouched. Limitation: the two
  windings of an `a` → `a` transformer cannot share one earth electrode; each side has its own (a shared electrode
  is listed under [Later work](#later-work)).
* **D15. `u1`/`u2` are winding voltages**, the rated voltages across the connection terminals: for example 22 kV
  for an `ab` primary, 12.7 kV for an `a` SWER secondary and 240 V for an `a` LV winding. Node `u_rated` stays
  $\sqrt3 \times$ the phase-to-earth voltage, so the nominal winding voltage implied by a node is $u_{rated}/\sqrt3$
  for a phase-to-earth and $u_{rated}$ for a phase-to-phase connection. Off-nominal ratios follow automatically and
  are legitimate, so C++ does not check `u1` against the node voltage; M9 implemented a warning-level Python rule
  ($|u/U_{nom} - 1| > 0.2$) for the classic mistake of `u1 = 22 kV` on an `a` winding at a 22 kV node, as a
  `ValidationWarning` (D42).
* **D16. General phase shift.** `phase_shift()` $= \varphi(\text{to}) - \varphi(\text{from}) + \text{clock}\cdot\pi/6$
  (see [Phase shift](#new-component-single_phase_transformer)). Output angles stay physical.
* **D17. Loading from the complex sum per side:** `loading` $= \max(|\sum_p s_{from,p}|, |\sum_p s_{to,p}|) / s_n$.
  The `transformer` rule ($\sum_p |s_p| / s_n$) would overstate the loading of a phase-to-phase winding by
  $2/\sqrt3$ even in normal operation.
* **D18. Working fixed taps.** `tap_side`, `tap_pos`, `tap_min`, `tap_max`, `tap_nom`, `tap_size`, `uk_min`,
  `uk_max`, `pk_min` and `pk_max` work as for `transformer` (the same tap-adjusted voltages and
  `tap_adjust_impedance`), and `tap_pos` can be updated through `TransformerUpdate`. Unlike `transformer`, all tap
  fields are optional: `tap_size` and `tap_nom` default to 0, `tap_pos` to `tap_nom` and `tap_side` to from, which
  gives the nominal ratio. Each tap limit is optional (M5b, so that a given `tap_pos` is not silently clamped to 0),
  with this rule:
  * The tap direction is decided once, at construction, from the provided limits. It is reversed (a higher tap
    position lowers the voltage of the tapped winding) when both limits are provided and `tap_max < tap_min`, when
    only `tap_min` is provided and `tap_min > tap_nom`, or when only `tap_max` is provided and `tap_max < tap_nom`.
    Otherwise it is normal, also without limits. For `tap_min == tap_max` this is normal, where `transformer` takes
    reversed; the tap position is then fixed at that value.
  * A provided limit always limits the tap position in that direction, at construction and in every update: in a
    normal tap changer `tap_min` is the lowest and `tap_max` the highest position, in a reversed one the other way
    round.
  * A limit that is not provided leaves its side unbounded, at construction and in every update. Not provided is
    remembered (stored as `na`); the `tap_min()`/`tap_max()` getters then return the current tap position $t$ when it
    lies on the side of `tap_nom` of that limit, and `tap_nom` otherwise (for example $\min(t, t_{nom})$ for `tap_min` of a normal tap
    changer and $\max(t, t_{nom})$ for a reversed one), so `tap_adjust_impedance` uses the provided limit on its own
    side and `uk`/`pk` on the other.
  * Without limits the tap position is not limited and the direction is normal.

  `uk_min`/`pk_min` only apply when `tap_min` is provided and `uk_max`/`pk_max` only when `tap_max` is provided;
  otherwise the short circuit parameters stay at `uk`/`pk` on that side of `tap_nom`. The tap regulator is listed under
  [Later work](#later-work); a `transformer_tap_regulator` on a `single_phase_transformer` raises `InvalidRegulatedObject`.

#### Taken in M6

The centre-tapped (split-phase) secondary of `single_phase_transformer` (all decided by Michael, 2026-10-16):

* **D19. One centre-tapped value, to side only.** `WindingConnection1P` gets only `ab_center_tapped = 7`, in C++ and
  Python. The legs map to slots a and b and the split-phase node is always a node of its own, so `ac`/`bc` variants
  would add test surface without a use; 8 and 9 stay unused (appending them later is ABI-compatible). The value is
  only allowed on the to side: on the from side it raises `InvalidPhaseConnection` with the ID ("a centre-tapped
  winding is only allowed on the to side"). 8, 9, 0 and `na` keep raising `MissingCaseForEnumError`. The primary is
  any of `a`, `b`, `c`, `ab`, `ac`, `bc`, as in M5.
* **D20. Half-winding impedances `uk_half` and `pk_half`** (optional). They are the short-circuit test between the
  primary and *one* half, with the other half open, on the pair rating $\min(s_n, s_n/2) = s_n/2$ and the half
  voltage $u_2/2$; the same convention as `three_winding_transformer`, whose `uk_12` is on
  $\min(s_{n,1}, s_{n,2})$. In per unit on $s_n$ and the half voltage: $r_h = 4\,p_{k,half}/s_n$,
  $|z_h| = 2\,u_{k,half}$ and $x_h = \operatorname{sgn}(u_{k,half})\sqrt{\max(|z_h|^2 - r_h^2, 0)}$.
  * **Defaults.** Not provided, each follows the split of Kersting (and Hopkinson), which OpenDSS's own centre-tap
    example uses: the primary $0.5R + j0.8X$ and each half $R + j0.4X$, in per unit on $s_n$, with $R + jX$ from
    `uk`/`pk`. That gives $z_{half} = 1.5R + j1.2X$, so `pk_half` defaults to $0.375\,p_k$ and `uk_half` to
    $\tfrac12|1.5R + j1.2X|$. The two default independently.
  * **On a winding that is not centre-tapped,** a non-empty `uk_half`/`pk_half` raises `InvalidPhaseConnection`
    ("… refers to a centre-tapped winding and must be left empty"), following D11: reject rather than ignore.
  * **Valid ranges** (implemented in M9 as Python rules, D41; no C++ check, as `transformer` does not check
    `uk = 0` either):
    `uk_half >= pk_half / (sn/2)`; `uk_half > uk/2` (otherwise the half branch $z_b$ is 0 or negative, and singular
    at equality); `pk/4 <= pk_half <= pk/2` (non-negative branch resistances).
* **D21. `u2`, `uk`, `pk` and `sn` describe the full winding.** `u2` is the full winding, end to end: 480 V for
  240-0-240 V, consistent with D15 (the voltage across the connection terminals a-b). `uk`, `pk` and `sn` are the
  full-winding nameplate: the primary against both halves in series, at rated `sn`. Each half is $u_2/2$. The
  split-phase node has `u_rated` $= 240\sqrt3 = 415.692$ V, because each leg is 240 V to the earthed centre tap; the
  nominal winding voltage implied by the node is therefore $2u_{rated}/\sqrt3$. M9 implemented this form in the
  warning rule of D15 (D42).
* **D22. Earthing of the centre tap.** Solid by default, or through `r_grounding_to`/`x_grounding_to`, the M5
  fields; not provided means solid. `r/x_grounding_from` follows the M5 rule for the primary. Caveat: PGM has no
  neutral node, so wye loads on the split-phase node are leg-to-*true-earth*. With an electrode, the whole
  out-of-balance load current therefore returns through the electrode. That is right for earth faults (the current
  to true earth returns through the electrode) and pessimistic for loads with a neutral conductor (MEN). Loads to a
  floating neutral are listed under [Later work](#later-work).
* **D23. Loading of a centre-tapped side:** $\max(|\sum_p s_{from,p}|,\ 2|s_{to,a}|,\ 2|s_{to,b}|)/s_n$. Each half
  is rated $s_n/2$: for example 12.5 kVA per leg on a 25 kVA unit is rated current at 120 V. Balanced full load gives
  1.0 on both sides, and all load on one leg at $s_n/2$ gives 1.0, where D17's $|s_a + s_b|$ would say 0.5. The
  primary keeps D17.
* **D24. Other points of the centre-tapped model.**
  * **Taps** are allowed, as in M5. A to-side tap scales both halves ($u_2'/2$). The half branch $z_b$ is fixed at
    the nominal tap and the primary branch takes the tap dependence of `uk`/`pk` (`uk_min/max`, `pk_min/max`):
    $z_0 = z_f(tap) - z_b/2$, because the tap winding is physically on the primary. Without `uk_min/max` this is the
    nominal split. There is no `uk_half_min/max`.
  * **Clock** 0 or 6. Clock 6 reverses both halves, so $U_a$ and $U_b$ change sign; the phase shift gets $+\pi$
    as in M5.
  * **Low susceptance:** never on a centre-tapped side (the centre tap is earthed). A phase-to-phase primary keeps
    it (D14).
  * **Status.** With the from side open, the primary is eliminated (the Schur complement in the 3×3 primitive): the
    halves stay coupled through $2z_b$, a **balancer**, which with $y_m = 0$ is
    $\frac{1}{2z_b}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$ in half-winding per unit, not zero. This differs from
    the "exactly 0" of a single connected M5 side, and is kept because it is physical. With the to side open, both
    halves are eliminated; this is exactly 0 when $y_m$ is below `numerical_tolerance` (the M5 rule), and otherwise
    the open-end shunt $\tfrac{y_m}{2} + (z_f + 2/y_m)^{-1}$ of M5.
  * **Magnetising:** the $\pi$ convention extended to $\mathrm{diag}(y_m/2, y_m/4, y_m/4)$, in per unit on $s_n$
    and the winding voltages. The two quarters at half voltage carry the same power as the secondary half of M5.
    OpenDSS places it differently, so the OpenDSS cases use $i_0 = p_0 = 0$.
  * **Short circuit at a split-phase node:** `single_phase_to_ground` on a or b (leg to earth, which is the centre
    tap when it is solidly earthed), `two_phase` ab (480 V) and `two_phase_to_ground` ab are valid. `three_phase`,
    and the default phases (bc for the two-phase faults), raise `InvalidShortCircuitPhases`.
  * **480 V loads** wait for the delta loads of M7.

#### Taken in M6, for start values

The start values of split-phase buses:

* **D25. Slot patterns for the start values** (decided by Michael, 2026-10-16). `Vector(T)` rotates the slots by
  0°/−120°/+120° and each bus has one scalar `phase_shift`, so a three-phase start puts slot b of a split-phase bus
  60° off ($a^2 U_a$ instead of $-U_a$).
  * **Type.** `enum class SlotPattern : std::uint8_t { three_phase = 0, split_phase = 1 }` in
    `<core>/calculation_parameters.hpp`. M7 reuses it for $U_{\Delta,nom}$: $\sqrt3$ for `three_phase`, 2 for
    `split_phase`.
  * **Incidence.** `BranchPhaseIncidence` gets `to_pattern` (default `three_phase`); only a
    `single_phase_transformer` with a centre-tapped to side reports `split_phase`. `all_phases()` is unchanged.
  * **Topology.** A step of `Topology::propagate_phases()` after the mask fixed point and its validation, so it only
    runs when `branch_phase_incidence` is non-empty. A worklist is seeded at the to node of every centre-tapped
    winding whose to side is connected and spreads a flag over the phase graph along the `per_phase` and `link` arcs
    only, never through a winding. The flag only grows, so the result does not depend on the visit order. A bus is
    `split_phase` iff it is flagged **and** has no slot c. A bus with slot c (for example an abc bus that is also
    fed three-phase) is neither flagged nor passes the flag on: otherwise the flag would leak back through the
    three-phase network into unrelated single-phase laterals (a refinement found by the topology tests during the
    implementation). The result is `MathModelTopology::bus_slot_pattern` (size `n_bus` iff `has_reduced_phase`,
    otherwise empty), read through `slot_pattern(bus)`.
  * **Solvers.** `detail::apply_slot_patterns<sym>(topo, u)` (`common_solver_functions.hpp`) sets
    $U_b = -U_a$ at every split-phase bus; slot a already carries the reference angle plus the phase shift of the
    bus, and slot c is left to the absent-slot mask that follows. It does nothing for a symmetric calculation or
    without `has_reduced_phase`. It is called in the Newton-Raphson flat and average-source starts and in the
    iterative-current start (always the average-source start). The Newton-Raphson linear start and the linear
    method solve the network directly and need nothing.
  * **Starts only.** The pattern changes start values, never the converged results, and the output angles stay
    physical. The one exception by construction is `linear_current`, which is a single iterative-current iteration
    and evaluates its loads at the start voltage: at a split-phase bus with `const_pq` or `const_i` loads it now uses
    the 180° start. Networks without a centre-tapped winding are unchanged.
  * **Phase shift.** $\varphi(\text{ab\_center\_tapped}) = 0$, so `phase_shift` $= -\varphi(\text{from}) +
    \text{clock}\cdot\pi/6$: 0 for `a` → ct and $-\pi/6$ for `ab` → ct (slot a at +30°, slot b at −150°).
  * **Limitation.** A transformer *fed from* a split-phase bus (for example `ab` → x or `b` → x off a 240-0-240 V
    node) gets the scalar phase shift of decision D16, which assumes the 120° pattern, so the start of its far side
    is 30° (`ab`) or 60° (`b`) off. This affects only the Newton-Raphson flat and average-source starts and the
    iterative-current start, and is rare: only loads and service lines normally hang off a split-phase bus. The
    default (linear) start of Newton-Raphson is not affected.
  * **Proof.** The native API test `API Model - centre-tapped single phase transformer` runs a loss-free `a` → ct
    network with `max_iterations = 1` and only a `const_pq` load without power (a network with only
    constant-impedance loads, or none, switches to the linear method). With the pattern, the flat and
    average-source starts are the exact solution and Newton-Raphson and iterative current converge in one
    iteration; without it each raises `IterationDiverge` (checked once by disabling the call). With unbalanced loads
    and without the pattern, Newton-Raphson from the flat and average-source starts even converged to a spurious
    solution (5 kV on a 240 V leg), which the same test also catches.
  * **Performance.** Gated as above: on a three-phase network the topology step does not run, `bus_slot_pattern`
    stays empty and the solver helper returns at once.

#### Taken in M7

Delta (phase-to-phase) loads (decided by recommendation under Michael's blanket approval, 2026-10-16):

* **D26. Enum and default.** `enum class LoadGenConnection : IntS { wye = 0, delta = 1 }` (C++ and Python),
  input-only. Not provided (`na`) means `wye`; any other value raises `MissingCaseForEnumError`, as for D10 and D12.
  * **Deviation from the design text:** the attribute also appears on `sym_load` and `sym_gen`. `LoadGenInput` is a
    single template input class, so the code generator emits one attribute list for the symmetric and asymmetric
    instantiations, and a separate asymmetric-only class would clash with the generated alias `AsymLoadGenInput`.
    `sym_load` and `sym_gen` accept `wye` and `na`; `delta` raises `InvalidPhaseConnection` with the ID ("a delta
    connection is only supported for asym_load and asym_gen"). The later delta `sym_load` reuses the attribute.
  * **C API compatibility.** Buffers made with `PGM_create_buffer` are uninitialised, so a caller that sets only the
    attributes it knows leaves `connection` as garbage: a value other than 0 or 1 fails construction with
    `MissingCaseForEnumError`, and a garbage value of exactly 1 would silently make an `asym_load` delta (softened by
    the experimental gate until the maintainers graduate the feature, D45). C API callers must therefore call `PGM_buffer_set_nan` on input buffers before
    setting attributes, as the C examples now do.
* **D27. Meaning of the entries.** For a delta load_gen, `p_specified`/`q_specified` are the **ab, bc, ca** values.
  Entry k is the branch from the first to the second letter: $U_{ab} = U_a - U_b$, $U_{bc} = U_b - U_c$,
  $U_{ca} = U_c - U_a$. A single phase-to-phase load is a delta load with one non-zero entry; a load between a and
  c uses the ca entry, also at an `ac` node. The orientation does not affect the power (flipping $U_k$ also flips
  $I_k$). Signs follow the existing load and generator reference directions.
* **D28. Absent phases (D3) for delta.** An entry is "specified" by the same exact test as for wye (not NaN and not
  0), and a specified entry needs **both** of its phases at the node: abc allows ab, bc and ca; ab (three-phase
  pattern or split-phase) allows ab; bc allows bc; ac allows ca; a, b or c allow none. Otherwise
  `LoadGenOnAbsentPhase` is raised, with the existing prefix (so the existing error mappings still match) and the
  tail "A delta (phase-to-phase) entry needs both of its phases: ab needs a and b, bc needs b and c, ca needs c and
  a." The wye message is unchanged. The check runs in the same place and under the same conditions as for wye, and
  the batch cache run still swallows it. As a defensive layer in the solvers, an entry is *active* iff it is not
  exactly 0 (NaN counts as active, as for wye) and both of its phases are present; inactive entries contribute
  exactly nothing.
* **D29. Voltage dependence and per unit.** The power base is the usual $S_{b,1p} = 10^6/3$ VA, so a delta entry
  is $S_k = \pm(p_k + jq_k)/S_{b,1p}$ through the existing `set_power`; the voltage base is the phase-to-earth base
  $u_{rated}/\sqrt3$ of the node, so $U_{\Delta,k} = (U_p - U_q)/U_b$. The nominal branch voltage comes from the
  slot pattern of D25: $U_{\Delta,nom} = \sqrt3$ (u_rated) for `three_phase` and $2$ ($2u_{rated}/\sqrt3$, 480 V
  on a 240-0-240 V node) for `split_phase`. The branch power is
  $S_k(U) = S_k (|U_{\Delta,k}|/U_{\Delta,nom})^n$ with $n = 0, 1, 2$ for `const_pq`, `const_i`, `const_y`. So
  `p_specified[k]` is the power at rated voltage across that branch (a `const_y` entry of 10 kW at a 400 V node is
  16 Ω). This matches the OpenDSS `kV` of a delta load, which is the branch voltage. Wye is the special case
  $D = I$, $U_{nom} = 1$.
* **D30. Output: per-phase node quantities.** `i` is the conductor current $|I_p|$ with $I = D^T I_\Delta$, and
  $s_p = U_p \overline{I_p}$ gives `p`, `q`, `s` and `pf` per phase. This fits `ApplianceOutput` unchanged, keeps
  the per-phase node power balance, keeps `calculate_source_result` and the voltage-regulator distribution (which
  sum `load_gen.i` or `.s` per phase) correct, and equals the OpenDSS terminal values (`CktElement.Powers`);
  $\sum_p p_p = \sum_k \mathrm{Re}\,S_k(U)$. The surprise is documented: at balanced voltages a resistive ab
  load P gives $p_a = p_b = P/2$, $q_a = -P/(2\sqrt3)$, $q_b = +P/(2\sqrt3)$ and 0 on c; at a split-phase node
  ($U_b = -U_a$) it gives $s_a = s_b = S/2$. Branch-level output (p_ab, ...) would need new solver output fields
  and can be added later.
* **D31. Short circuit.** Loads are ignored in short circuits (`get_null_sc_output`, and the short-circuit solver
  has no load_gen code), so delta loads need no short-circuit change. A model whose only special feature is delta
  loads may still use the **symmetric** short circuit (all faults three-phase). The experimental gate still applies.
* **D32. Where delta is allowed: at all nodes, abc included.** `ComponentTopology::load_gen_connection` and
  `MathModelTopology::load_gen_connection` are **empty unless a delta load_gen exists** (per math model, and only
  for a connected delta load_gen). The solvers test `topo.is_delta_load_gen(i)`, which is false when the vector is
  empty, and the wye arithmetic is untouched, so wye-only networks keep the bit-identical fast path.
* **D33. Calculation gates.**
  * **(a) Experimental** until the maintainers graduate the feature (D45): any delta load_gen in the input (status and connectivity ignored, like M4)
    raises `ExperimentalFeature` ("A delta-connected load or generator is an experimental feature"). This check
    comes after the reduced-phase gate, whose message wins.
  * **(b) State estimation** raises `ReducedPhaseNotSupported` for both symmetries and all methods: state
    estimation is listed under [Later work](#later-work), and the measured-values load/gen distribution has no delta model.
  * **(c) Symmetric power flow** raises `ReducedPhaseNotSupported` (open question 2 stays out of scope). The
    symmetric short circuit is allowed (D31).
  * **(d) One exception class.** Both use `ReducedPhaseNotSupported` (D6) with a second message flavour
    (`ReducedPhaseNotSupported::Cause::delta_load_gens`): "{feature} is not supported in a network with
    delta-connected loads or generators." The C++ validation regex accepts both flavours (the Python mapping
    follows with the Python enum).
  * **(e) Voltage regulators.** A `voltage_regulator` on a delta load_gen raises the existing
    `UnsupportedLoadGenTypeForVoltageRegulator`. A non-regulated delta load at a PV bus is allowed.

#### Taken in M8

Short-circuit phase validity (decided by recommendation under Michael's blanket approval, 2026-10-16):

* **D34. Model-level check on the fault phases.** `detail::check_fault_phases_present(state)` in
  `<core>/main_core/calculation_input_preparation.hpp` runs in `prepare_short_circuit_input` when a math model has
  `has_reduced_phase`, before the shunt check. For every fault with status on whose node is in a math model with
  `has_reduced_phase`, the resolved fault phases (`get_fault_phase()`, mapped by `fault_phase_mask` in
  `<core>/calculation_parameters.hpp`) must all be present at the node (`topo.bus_phase`, after topology and phase
  propagation). A fault with status off (as in D3) and a fault at an isolated or de-energised node are skipped; the
  latter keeps the null output (i_f = 0), as in three-phase networks. The first failing fault in input order raises.
  The M2 solver check stays unchanged as an ID-less backstop for direct solver users; through the model it can no
  longer be reached.
* **D35. Reuse `InvalidShortCircuitPhases` with a new constructor.**
  `InvalidShortCircuitPhases(fault_id, node_id, fault_phases, node_phases, default_phases)` gives "The short circuit
  phases (b) of fault 601 are not all present at its node 4 (node phases: a); this short circuit type cannot be
  applied there", with the tail "; the phases are the default of the fault type, so set fault_phase explicitly" when
  `fault_phase` was not given or is `default_value` (`Fault::has_default_phase()`) and some explicit phase of that
  fault type fits the node: `single_phase_to_ground` needs at least one node phase, `two_phase` and
  `two_phase_to_ground` at least two, `three_phase` never (it always needs `abc`, so a hint would mislead). No new class and no Python
  change: the message keeps "short circuit type", so the Python mapping and the C++ validation regex still map it
  to `InvalidShortCircuitPhaseOrType`, as in M2 and M6. The two existing messages are unchanged.
* **D36. Default phases are not adapted to the node.** `default_value` keeps resolving statically
  (`single_phase_to_ground` → a, `two_phase` and `two_phase_to_ground` → bc, `three_phase` → abc). The default
  single-phase fault is valid at an `a` SWER node and raises, with the default tail, at a `b`- or `c`-only node; at
  a split-phase node the default single-phase fault is valid and the default two-phase faults raise (as M6
  documents). Not resolving per node, because node phases follow status updates (the faulted phase would change
  silently between scenarios), because the one-type-and-phase rule compares resolved phases (two default faults at
  an `a` and a `b` node would fail confusingly), and because rejecting first is forward-compatible. Users should
  use `to_connection = a` for SWER secondaries so that the default applies.
* **D37. Batch: the cache run swallows `InvalidShortCircuitPhases`.** `JobAdapter::cache_calculate_impl` catches it,
  consistent with D3 and D4: `fault_phase`, `fault_object` and `status` can be updated, and topology changes per
  scenario move the node phases, so the update data may fix the base fault. Each scenario runs the check on its own
  and fails on its own, with the fault ID. During a calculation the class is only thrown by D34 and the solver
  backstop (`Fault::check_sanity` throws it at construction and update, which the cache run does not do).
  `InvalidShortCircuitPhaseOrType` for mixed types in the base input keeps its existing, unswallowed behaviour.
  Note: when the base input has mixed fault types and also a fault on an absent phase, the model check
  (`InvalidShortCircuitPhaseOrType` for the mixed types) runs first and the cache run swallows the phase error. Each
  scenario then fails on its own with `InvalidShortCircuitPhaseOrType`; with mixed types only, the whole batch still
  fails at once.

#### Taken in M9

Python validation and documentation (decided by recommendation under Michael's blanket approval, 2026-10-17):

* **D38. What Python checks and what stays in C++.**
  * **Python checks** the rules on the data of each component, plus simple ID lookups (the node `u_rated` behind a
    winding, the load behind a voltage regulator).
  * **C++ keeps** every rule that needs the derived node phases or a calculation: D1 (`MissingPhaseForWinding`,
    `DanglingConductor`), D3 and D28 (`LoadGenOnAbsentPhase`, which includes the delta entries against the node
    phases), D4 (`ShuntAtReducedPhaseNode`), D34 and D36 (fault phases against the node phases) and the calculation
    gates (symmetric calculations, state estimation, `ExperimentalFeature`).
  * **Why.** The node phases come from propagation from the sources through lines, links and windings, which depends
    on status and has the "all terminals of a winding" rule. Copying it into Python would give two implementations
    that can drift apart. The C++ errors already name the ID and run once per scenario. The Python framework has no
    topology checks today (not even `ConflictVoltage` or connectivity), and a partial check would give false
    confidence. The expected-error cases `delta-on-absent-phase`, `swer-fault-on-absent-phase` and
    `missing-phase-for-winding` are pinned by a test: they pass the Python validation and still fail in the
    calculation. A best-effort Python check of the node phases is listed under [Later work](#later-work).
* **D39. `asym_line.phases` follows D11.**
  * The present conductors are the phases in `phases`, plus the neutral `n` when `r_nn` or `x_nn` is given. Every
    matrix entry whose two conductors are both present is required (`MissingValueError` per field), which covers r, x,
    the c matrix and the neutral entries.
  * A value on an absent phase, on a missing neutral, or `c0`/`c1` on a one- or two-phase line raises the new
    `NotApplicableValueError(MultiFieldValidationError)`, with the fields `[phases, <field>]`.
  * An unknown `phases` value raises `InvalidEnumValueError`, and the row gets no per-phase checks.
  * The upstream abc rules (the neutral subset, the c-matrix subset, `c0`/`c1`, "not all c missing") now run on rows
    where `phases` is not provided or is `abc`. The rules "self value > 0" and "mutual r/x >= 0" ignore NaN and run
    on every row.
* **D40. Drop the sign rule on the mutual capacitance** (a pre-existing defect, see
  [above](#pre-existing-defects-found-during-this-investigation)). C++ uses the c matrix as the nodal (Maxwell) matrix,
  so any sign is valid; positive values, as in the upstream `rmatrix-xmatrix-cmatrix` data, stay valid. Without
  this, 25 of the 33 reduced-phase cases with a single-phase transformer, `phases` or a load `connection` fail
  the validation on this rule alone. It is a commit of its own, so that it can be proposed upstream separately, like M0.
  The self capacitance must still be positive.
* **D41. `validate_single_phase_transformer`.**
  * **Required:** the branch fields plus `u1`, `u2`, `sn`, `uk`, `pk`, `i0`, `p0`, `from_connection` and
    `to_connection`. `clock`, the tap fields, `uk_half`/`pk_half` and the earthing fields are optional (D13, D18, D20,
    D14).
  * **Ranges**, as for `transformer`: `u1`, `u2`, `sn` > 0; `0 < uk < 1` and `uk >= pk/sn`; `pk >= 0`;
    `p0/sn <= i0 < 1`; `p0 >= 0`; `uk_min/max` and `pk_min/max` with the same rules and defaults.
  * **Connections** must be valid `WindingConnection1P` values, and `ab_center_tapped` on the from side gives
    `InvalidValueError` (D19).
  * **Clock** other than 0 or 6 modulo 12 raises the new `SinglePhaseTransformerClockError`; the existing
    `TransformerClockError` talks about wye windings and does not fit (D13).
  * **Earthing:** `r/x_grounding_<side>` on a phase-to-phase winding gives `NotApplicableValueError`, and
    `r_grounding_*` must be $\geq 0$ (D14).
  * **Half windings:** `uk_half`/`pk_half` on a winding that is not centre-tapped gives `NotApplicableValueError`.
    On a centre-tapped winding, with the C++ defaults applied to whichever value is not given: `uk_half > uk/2`,
    `pk/4 <= pk_half <= pk/2` and `uk_half >= pk_half/(sn/2)` (D20).
  * **Taps (D18):** `tap_side` must be a valid `BranchSide` and `tap_size >= 0`. The direction comes from the limits
    that are given. `tap_pos` (default `tap_nom`) must lie within the given limits, each in its own direction, and an
    omitted limit is unbounded; `tap_nom` must lie between the two limits when both are given. `tap_pos` can be
    updated, so `validate_batch_data` checks it per scenario. This mirrors `transformer`, where C++ silently clamps
    the value instead.
* **D42. The u1/u2 rule is a warning, not an error** (D15, D21).
  * The framework had no warning level. `class ValidationWarning(UserWarning)` is added to
    `power_grid_model.validation` and raised with `warnings.warn`.
  * **Trigger:** $|u/U_{nom} - 1| > 0.2$, where $U_{nom}$ is $u_{rated}/\sqrt3$ for `a`, `b`, `c`; $u_{rated}$ for
    `ab`, `ac`, `bc`; and $2u_{rated}/\sqrt3$ for `ab_center_tapped`.
  * **Once per call** of `validate_input_data` or `validate_batch_data`, on the input data: `u1`, `u2` and `u_rated`
    are input-only, so a per-scenario loop would only repeat it.
  * **Why not an error:** off-nominal ratios are legitimate (D15), so an error would block valid models. The warning
    does not change the return value or `assert_valid_*`.
* **D43. Load `connection`.** The value must be a valid `LoadGenConnection` on all four load_gen types
  (`InvalidEnumValueError`); `delta` on `sym_load`/`sym_gen` gives `InvalidValueError` (D26); a `voltage_regulator` on
  a delta `asym_load`/`asym_gen` gives `InvalidVoltageRegulationError`, the error the const_power rule already uses
  (D33e), at most one per regulator.
* **D44. Sensors can measure a `single_phase_transformer`.** It is added to the branch lists of the power and current
  sensors (the general list, `branch_from` and `branch_to`). C++ accepts any `Branch`, a power flow ignores sensors,
  and state estimation raises `ReducedPhaseNotSupported` anyway. A tap regulator on a `single_phase_transformer`
  already gives `InvalidIdError`, because `validate_regulator` only lists `transformer` and
  `three_winding_transformer`; a test pins this.
* **D45. The experimental gate stays (D9, D33a).** Only the wording that tied the gate to M9 changed, to "until the maintainers
  graduate the feature". The pull request proposes removing it.
  * Graduating a feature is the decision of the maintainers (`terminology.md`: undefined behaviour, no guarantees).
  * Real limitations remain: state estimation, symmetric calculations, sources and voltage regulators on
    reduced-phase buses; the start-value limitation of D25; the iterative-current divergence on a floating common
    mode with only const_pq or const_i loads (const_y loads anchor it since M11, D61); and the topology rules that
    only C++ checks (D38).
  * **The public Python API has no experimental switch:** `calculate_*` does not accept `experimental_features`, only
    the private `_calculate_*` methods do (the tests use them). Removing the gate is the only way to make the
    feature publicly usable from Python, which is an API decision for upstream.
  * Removing the gate is cheap: one check in `check_no_experimental_features_used` in
    `<core>/main_model_impl.hpp`, plus the documentation. The user manual states how to enable it today: the C API
    option `PGM_set_experimental_features`, or the private `_calculate_power_flow(...,
    experimental_features="enabled")` in Python, flagged as private and subject to change.
* **D46. Documentation: a new user-manual page with a worked Python example that runs in a test; no notebook**
  (Michael's earlier call). `docs/user_manual/reduced-phase-networks.md` replaces the "example notebook" of the M9
  row, and `tests/unit/test_docs_reduced_phase_example.py` executes its code blocks.

#### Taken in M10

Validation against two real networks (decided by recommendation under Michael's blanket approval, 2026-10-18; D47 was
reported to him, because it changes the approved plan):

* **D47. Maximum demand is 20 % of the transformer kVA, minimum demand 10 %, power factor 0.95.** The 60 % first
  approved does not work: OpenDSS does not converge on either network and the voltages collapse (SMR8 down to 0.29
  pu, KLO14 to 0.53 pu), because 60 % of the installed transformers (40 MVA and 35 MVA) is 24 to 26 MW on rural 22 kV
  feeders, and the regulators would run at 153 % to 281 % of their rating. At 20 % every regulator is within its
  rating, the feeder peaks are 8.1 MW (SMR8) and 7.0 MW (KLO14), and the lowest node voltage is 0.845 pu (SMR8, on a SWER
  node) and 0.910 pu (KLO14, on a split-phase LV busbar; its lowest SWER node is at 0.933 pu). The generator has a
  voltage gate of 0.8 pu and stops with a message that names this decision when the loads are too high. `MAX_DEMAND` and `MIN_DEMAND` are constants in the generator.
* **D48. The ISO coincidence factor applies at maximum demand only.** The SWER loads behind an isolating transformer
  (ISO) get $\min(1, S_{ISO} / (0.2 \sum S_{tx}))$, with $S_{tx}$ the rating of the SWER distribution transformers behind it. The
  factor models the diversity of the peak. Minimum demand is base load, and every ISO is within its rating without it. At 20 %
  it binds on two ISOs of SMR8 (NISSEN_ISO 0.625, BUNGARRA_ISO 0.9346) and on none of KLO14.
* **D49. Voltage regulators are fixed-neutral `YNyn0` banks.** `backend.py` builds each regulator as three
  autotransformers under a `RegControl` that starts at tap 0. At the neutral tap this is a 1:1 connection plus the
  leakage impedance, so both tools get a `transformer` (clock 0, grounded wye on both sides) with the rating, `xhl` and
  `%loadloss` of the sheet. Regulator action (tap changing) is not part of the validation. This includes
  PLENTY_VLY_RISE of KLO14 (22/22 kV, 200 kVA, without `_REG` in its name).
* **D50. ISOs are two-winding.** The ISO is a `single_phase_transformer` from a phase pair of the 22 kV side
  (`Conn_Type`; clock 6 if the pair is listed in reverse order) to phase a of the SWER network. The open third winding
  of `backend.py` is dropped. The OpenDSS buses are written literally from `Conn_Type`, so OpenDSS checks the clock-6
  mapping independently.
* **D51. Split-phase transformers use an explicit Kersting split.** The 250-0-250 V transformers are centre-tapped
  (`ab_center_tapped`, $u_2$ = 500 V) with the Kersting defaults of D20. OpenDSS gets the matching `XHT` and `XLT`
  instead of the defaults of `backend.py` (`XHT` = 35 %, `XLT` = 30 %), which are not the Kersting split. The YPrim of
  every distinct parameter set is checked against the PGM model.
* **D52. The 44 phase-to-earth units of KLO14 are 12.7 kV / 250 V.** The data gives them a 22/0.433 kV rating but wires them from phase R
  to earth. Both ratings are taken phase-to-earth, which keeps the turns ratio of the data. They are
  `single_phase_transformer`s, like the SWER units.
* **D53. $i_0 = p_0 = 0$ everywhere.** PGM places the magnetising branch as a $\pi$ (half on each side), OpenDSS lumps it
  at winding 1, as in M5. `uk` and `pk` come from the sheet's `xhl` and loadloss, and OpenDSS gets the matching
  `XHL` and `%loadloss`. The low-susceptance zero-sequence regularisation of the delta side of a `transformer`
  (and of a phase-to-phase `single_phase_transformer` winding, D14) is mirrored by shunt `Reactor`s, whose power is
  added to the from side of the reference, as in M5 (the reactors equal PGM's value to 1.0e-6, the `eps` regularisation).
* **D54. Near-zero connectors are merged, and the zone-substation $u_k$ is raised to 1e-4.** The data has connector
  lines (LineCode 0, 0.001 ohm/km, and others) with a series impedance down to 1.4e-7 ohm, and a zone substation with
  $u_k$ = 1.41e-6. They make the admittance matrix so ill-conditioned that the PGM Newton-Raphson stalls at a maximum
  deviation of about 3e-7, above its error tolerance of 1e-8, also after 1000 iterations; PGM and OpenDSS agree to only about 2e-7.
  * **This is pre-existing, not caused by the reduced-phase code.** The model with only its three-phase parts (every
    `single_phase_transformer`, `asym_line` and their nodes and loads removed; the upstream code path) stalls at the
    same level (0.8e-7 to 5.3e-7).
  * **A floor on the impedance was tried and rejected.** M10b lengthened every line below 1e-2 ohm. That worked for PGM,
    but it changed the OpenDSS solution of the source network by up to 0.4 % in voltage magnitude (5.3e-3 and 6.6e-3 as phasors) and by +4 % to +10 % in
    losses, so the cases no longer described the real networks.
  * **Decision.** Every line with $|z_1| L$ below 1e-3 ohm is removed and its end nodes are merged (union-find; the
    representative has the fewest branches to the source, then the lowest id), in the input of both tools. Lines,
    transformer terminals and capacitors move to the representative. The generator asserts that a merge never joins a 22 kV node and a SWER node. SMR8 loses 887 of
    its lines (4329 to 3442 nodes), KLO14 1310 (4533 to 3224 nodes); the smallest line left is 1.01e-3 ohm.
  * **The zone-substation $u_k$** is raised to 1e-4 (0.00052 ohm on the 22.9 kV side), because after the merge it is
    the smallest impedance left: with 1e-5 or 3e-5, KLO14 misses the tolerance rule (ratio 0.145, on the power through the
    200 kVA regulator).
  * **Deviation from the source network** (every line, the sheet $u_k$, solved by OpenDSS with the same tolerance): at most 5.5e-5 in node voltage
    (phasor metric, merged nodes at their representative) and at most 0.1 % in losses; the feeder power changes by at most 2.3e-5. The README of each case lists the
    change per node class.
  * Faults sit on kept nodes only, so that the fault node exists in both tools.
* **D55. The references carry no three-phase `line`, short-circuit `node` or `asym_line` output.** Whole component
  types can be left out of a reference (the harnesses compare the keys of the reference), but not a subset of rows. The
  power-flow reference has `node`, `asym_line`, `single_phase_transformer` and `transformer`; the short-circuit
  reference has `fault`, `single_phase_transformer` and `transformer`. The node voltages at both ends pin every line, and the
  size budget (about 10 MB for the new data) did not allow more: `asym_line` output would add 1.05 MB to the short-circuit cases and
  `node` output about 5 MB. The datasets use the compact JSON format (a root `attributes` list, one list row per line), and references
  are rounded to 12 significant digits.
* **D56. Tolerances.** The rule is that of the earlier cases: the smallest round values with at least 10x margin over all
  scenarios, with a `dict` for `atol` where one field needs its own value. The power `atol` is capped at 1e-3 W (var),
  so the power flow gets rtol 1e-6; rtol 1e-7 would need 1e-2 W.
  * **Power flow:** rtol 1e-6, atol 1e-8, and 1e-3 on `[pq]_(from|to)`. Node voltages agree to 1.1e-10 (SMR8) and 9.8e-11 (KLO14)
    relative, the branch powers to 1.6e-3 W and 1.2e-3 W absolute. A branch power is computed from the voltages at its ends, so a branch
    with a small impedance (down to 1e-3 ohm) amplifies their round-off to an absolute error that no relative tolerance covers; the
    dict `atol` is limited to the power fields. The worst ratio of measured deviation to tolerance is 0.066 and 0.069. `newton_raphson` runs on both networks (3 and 4
    iterations), `iterative_current` on SMR8 only (8 and 12 iterations); on KLO14 it misses the rule (ratio 0.207, on the power of
    the regulator PLENTY_VLY_RISE). *M11 update:* both methods run on both networks, with iterative current at an error
    tolerance of 1e-10. SMR8 tightens to rtol 1e-7 (worst ratio 0.050 for Newton-Raphson, 0.044 for iterative current);
    KLO14 keeps rtol 1e-6 and iterative current gets its own `atol` of 1e-2 on `[pq]_(from|to)` (D62).
  * **Short circuit:** rtol 1e-8, atol 1e-8, and 1e-6 on `i_(from|to)` and its angle. Fault current magnitudes agree to 6.7e-11 and 9.2e-11 relative, and
    branch currents above 1 A to 6.0e-9 relative. Currents up to 10 A agree to 3.6e-8 A absolute: they are small currents through low-impedance branches,
    driven by small voltage differences, which have the same round-off amplification. rtol 1e-9 is not enough (the fault angles agree to 3e-10 only). The worst ratio is 0.033 and 0.030
    (the 66 kV current of the zone substation for a remote SWER fault).
* **D57. Licence sidecars are `MPL-2.0 AND BSD-3-Clause`.** The data of the cases derives from the Team-Nando data (Copyright (c) 2024,
  Team Nando - Public Repositories, BSD 3-Clause), so every data file of a case has a `.license` sidecar, and each `README.md` a header, with both
  copyright lines and the licence expression `MPL-2.0 AND BSD-3-Clause`. `LICENSES/BSD-3-Clause.txt` is added for `reuse lint`. There is no
  `REUSE.toml` override. The workbooks are not part of the repository: the generator downloads them from the pinned commit (or reads
  `--source-dir`) and checks their sha256.
* **D58. The OpenDSS reference follows what OpenDSS solves.**
  * **One-phase lines use $Z_1$ and $C_1$ only.** For a `LineCode` with `nphases=1`, a `Line` solves with the series impedance
    $r_1 + jx_1$ and the shunt $C_1 = b_1/\omega$; $r_0$, $x_0$ and $b_0$ are ignored. The matrix that the `LineCode` object reports is
    $(2Z_1 + Z_0)/3$ and must not be read. The `asym_line` of a SWER section with phase `a` therefore gets $r_{aa} = r_1 L$,
    $x_{aa} = x_1 L$ and $c_{aa} = b_1 10^{-6} L/\omega$ ($b_1$ in µS/km), and has no earth-return reactance, as in `backend.py`.
  * **Tolerance 1e-11, `MaxIterations` 200.** `Set Tolerance=1e-12` stalls at round-off on SMR8 and gives the same voltages as 1e-11. The earlier cases use
    1e-12; the README of each case states the change.
  * **Taps are written explicitly** (`taps=[t 1]`), because the `mintap`/`maxtap` of `backend.py` do not clamp. The unreduced source network is also solved, for the
    README table of D54.
  * **The generator checks before it writes anything:** the LineCodes against the PGM line rows (the $(2Z_1 + Z_0)/3$
    form), the YPrim of every distinct single- and three-phase transformer set against the numpy PGM model, the mirror reactors, the capacitors, the convergence, that the loads draw their kW + jkvar,
    the 0.8 pu voltage gate (D47) and the topology and counts.
* **D59. Runtime cost in the debug test run.** In the debug `ctest` (clang with the sanitizers), `Validation test batch - power flow` goes from 98.8 s to
  235 s (the two cases take 79 s and 58 s, of which `iterative_current` on SMR8 is 19 s), and `Validation test batch - short circuit` from 24.8 s to
  70.9 s. The time is spent in the sanitizer build on networks of 3.2k to 3.4k nodes. The merge of D54 cut the power-flow
  time from 256.5 s, which it took with the impedance floor. *M11 update:* the power-flow batch takes 238.6 s with
  iterative current on every M3 to M9 case and 256.6 s with iterative current on KLO14 as well (about 18 s more; debug
  build, one run each).

#### Taken in M11

The iterative current method (IC) runs on every power-flow validation case (decided by recommendation under Michael's
blanket approval, 2026-10-18; D62 is the one place where a tolerance is relaxed for one method, and Michael is asked to review it):

* **D60. Per-method solver options in the validation parameters.** `extra_params.<method>` may hold `error_tolerance`,
  `max_iterations`, `rtol` and `atol`. The C++ harness (`construct_case`) and `reduced_phase_compare_pgm.py` now read
  them, as the Python harness (`tests/unit/utils.py`) already did. A key in `extra_params.<method>` replaces the key of `params.json` as a whole, so
  a nested object such as a dictionary `atol` is replaced, not merged; the C++ harness did a recursive merge before
  and now matches Python. The defaults are `error_tolerance` 1e-8 and `max_iterations` 20. An error tolerance for all
  methods in one case is ruled out, because Newton-Raphson stalls at about 5e-11 on the SMR8 and KLO14 networks and
  does not converge at 1e-11 or below.
* **D61. The IC changes are gated on "some bus has an absent phase".** The gate is an asymmetric calculation where
  `topo.bus_phase(bus)` is not full for at least one bus; `has_reduced_phase` alone is not enough. Three-phase networks,
  abc networks with delta loads and a network flagged reduced-phase with all buses abc keep the old path, bit for bit.
  The gate is per math model, that is per island with its own source: an island whose buses are all abc keeps the old
  path even when another island of the same network has an absent phase (`has_reduced_phase` is set per math model).
  The M11 review checked this bit for bit with the abc delta case as one island next to the `phase-to-phase`,
  `split-phase` and delta `two-phase` cases. In a gated network:
  * **const_y loads are in the matrix.** The wye and delta const_y `load_gen`s are added to the diagonal of the
    prefactorised matrix as constant admittances and skipped on the right-hand side. A const_y load to earth then
    anchors a floating common mode (the `ab` secondary of an `ab` → `ab` transformer), where the explicit form
    diverged with a gain of about 1e6 per iteration. The per-bus const_y admittance of the input is cached and the
    matrix is refactorised when it changes (exact comparison), so a batch does not depend on the scenario order. The
    M10 networks have no const_y load and never refactorise. No validation batch changes const_y loads at an absent
    phase, so the M11 review checked it in Python: batches that change the const_y loads of seven reduced-phase and
    delta cases give bit-identical results with threading -1, 0, 1 and 2 and equal single calculations bit for bit
    (IC, `linear_current` and Newton-Raphson), and so do the SMR8 and KLO14 batches.
  * **`linear_current` treats const_y exactly.** It is the first iteration of IC, so the const_y loads are in its matrix as well.
    No validation case notices: the `-linear` cases are all const_y and run the `linear` method. Before this, `phase-to-phase`
    `linear_current` returned 7e5 pu garbage. The unit test "linear current at a network with an absent phase" checks
    that it equals the `linear` method when only const_y loads carry power.
  * **NaN of any bus and phase reaches the guard of D64.** `fold_max_deviation` replaces `std::max` and `max_val`
    in the Newton-Raphson and IC `iterate_unknown`, which could drop a NaN depending on the argument order.
    Finite deviations are unchanged.
  * **Measured bit-identity.** A capture of 21506 result keys (144 cases, 74 batch runs) before and after: every
    Newton-Raphson, linear, symmetric and three-phase key is identical. Only the IC and `linear_current` keys of absent-phase cases with const_y
    loads differ (IC node voltages by up to 2.5e-10 relative), and `phase-to-phase` IC no longer raises.
  * **Applying it to all networks** is listed under [Later work](#later-work), because it changes three-phase results.
* **D62. KLO14 runs IC with an IC-only `atol` of 1e-2 on `[pq]_(from|to)`.** *Taken under blanket approval; Michael
  is to review it, because it is a looser tolerance than Newton-Raphson has.* Newton-Raphson keeps 1e-3.
  * **What fails.** At the default IC the worst entry is the 795 W phase-c flow through the 200 kVA regulator PLENTY_VLY_RISE
    (ratio 0.207 with the 1e-3 `atol`). Raising `error_tolerance` to 1e-10 lowers it to 0.111, no further: the fixed
    point of IC solves the rounded matrix, and the entry is off by about 2e-4 W of round-off (Newton-Raphson 1.2e-4 W).
    The failure is round-off at the 10x margin (0.111 against 0.1), not convergence.
  * **Result.** With `error_tolerance` 1e-10 and the 1e-2 `atol` the worst ratio is 0.018 and 0.016 (9 and 12 iterations);
    node voltages and the other fields keep the tolerances of Newton-Raphson. The other option is to leave KLO14 Newton-Raphson only.
  * **C1 was evaluated first and rejected.** C1 is the residual (defect-correction) form of IC: iterate
    $\Delta U = M^{-1} (I(U) - Y U)$, $U \leftarrow U + \Delta U$, with $M$ the prefactorised matrix and $Y$ the true admittance
    matrix, so that the fixed point is $I(U) = Y U$ evaluated with the true matrix. It was implemented behind the absent-phase gate (kept as a patch in the
    local notes, not in the repository) and measured on KLO14 and SMR8:
    * It did not give a better accuracy. Against OpenDSS KLO14 reached 0.085/0.063 at error tolerance 1e-10 and
      0.058/0.050 at 1e-9 (passing the rule), but IC against Newton-Raphson was no closer: 0.049/0.030 (1e-10) and
      0.040/0.065 (1e-9) with C1, against 0.044/0.046 (1e-10) and 0.041/0.046 (1e-14) without it.
      C1 does not remove a systematic offset; it makes IC limited by the residual round-off, like Newton-Raphson.
    * It added a stall floor. With C1, IC stalls at a deviation between 3e-11 and 5e-11 on KLO14 and between 2e-11 and 3e-11 on SMR8
      (error tolerance 1e-11 raises), so 1e-10 leaves a margin of 2x to 3x only, where the form without C1 converges
      down to 1e-14.
    * The ratios move with the error tolerance (KLO14 scenario 0: 0.133, 0.058, 0.085 at 1e-8, 1e-9, 1e-10), which is noise.
      Other cases got worse: SMR8 at 1e-10 from 0.0285 to 0.0475, `swer-chain` at 1e-8/1e-8 from 3.6 to 10.3.
    * A small-grid check against a long-double reference (a chain with a 1e6x branch) showed the same: the form
      without C1, with C1 and Newton-Raphson all sit at the floor of machine epsilon times the condition number,
      so no unit test could separate C1 from the form without it.
* **D63. The IC cases run with an IC-only `error_tolerance` of 1e-10, and the tolerances are tightened.** The constant
  `IC_ERROR_TOLERANCE` of `reduced_phase_opendss_cases.py` is passed through `extra_params.iterative_current` in every case that lists IC, and
  `australian_swer_networks.py` imports it. At the default 1e-8 IC stops on its convergence residual (worst ratio 2.645 on
  `split-phase-flat-start`); at 1e-10 it agrees with OpenDSS as closely as Newton-Raphson (0.034). It is 9 to 14 iterations,
  below the limit of 20. Tolerances were then re-measured by the rule of the earlier cases (the smallest round `rtol` with `atol` ≤ 1e-6, then the smallest `atol`, worst ratio < 0.1 over all methods):

  | Case | Before (rtol / atol) | After | Worst ratio |
  |------|----------------------|-------|-------------|
  | `asym-line/single-phase-a` | 1e-6 / 1e-8 | 1e-8 / 1e-8 | 0.034 |
  | `asym-line/two-phase-bc-neutral-opendss-kron` | 1e-7 / 1e-8 | 1e-8 / 1e-8 | 0.034 |
  | `asym-line/swer-a` | 1e-7 / 1e-7 | 1e-8 / 1e-7 | 0.039 |
  | `asym-line/reduced-phase-feeder` | 1e-6 / 1e-6 | 1e-8 / 1e-6 | 0.080 |
  | `single-phase-transformer/swer-chain` | 1e-7 / 1e-5 | 1e-8 / 1e-6 | 0.049 |
  | `single-phase-transformer/split-phase` | 1e-6 / 1e-8 | 1e-8 / 1e-8 | 0.034 |
  | `asym-load-delta/split-phase` | 1e-6 / 1e-8 | 1e-8 / 1e-8 | 0.034 |
  | `asym-load-delta/two-phase` | 1e-6 / 1e-8 | 1e-7 / 1e-8 | 0.033 |
  | `australian-swer/smr8` | 1e-6 / 1e-8 (1e-3 on the branch powers) | 1e-7 / 1e-8 (same) | 0.050 |

  Every tolerance that is still looser than 1e-8 / 1e-8 is needed by Newton-Raphson too (near-zero line powers, nearly cancelling
  line reactive power). The `*-linear` cases are unchanged. Only `params.json` and the READMEs changed; a second
  generation is byte-identical.
* **D64. A NaN deviation raises `IterationDiverge`.** The loop of Newton-Raphson and IC runs while `max_dev > err_tol`; a
  NaN fails that comparison, so the loop ended as if it had converged and returned NaN results without an error.
  `IterativePFSolver::run_power_flow` now throws `IterationDiverge` as soon as the deviation is NaN. Finite results are
  unchanged (the 21506-key capture is identical). The unit test uses a three-phase bus whose zero sequence is held to
  earth by 1e-8 only, with an explicit const_y load: the voltages overflow and turn NaN within 100 iterations. A floating
  mode with only const_pq loads never gives NaN: IC stays finite (deviation about 5e5) and raises at `max_iterations`.
* **D65. IC keeps ignoring `calculation_initialization`.** IC always starts from the average-source voltage (with the slot patterns
  of D25); the option applies to Newton-Raphson only. The two `split-phase-flat-start` cases run IC from the
  average-source start, which equals the flat start there (8 and 9 iterations). Making IC honour the option is listed under [Later work](#later-work).
* **D66. Scope of IC in the validation cases.** Every reduced-phase, split-phase, delta and real-network power-flow case that lists `newton_raphson` also lists `iterative_current`:
  `phase-to-phase`, both `split-phase-flat-start` cases, KLO14 and the four expected-error cases (they raise the same exception type and message as
  Newton-Raphson). The `-linear` cases do not, because with all loads const_y PGM runs the `linear` method for every method name.
  The `*-linear` cases list `linear_current`, whose const_y treatment in absent-phase networks changed (D61).

### Decision index

| Decision | Subject | Section |
|----------|---------|---------|
| D1 | Winding without terminals at a connected side is an error | [Taken in M2](#taken-in-m2) |
| D2 | `PhaseCoupling::link` | [Taken in M2](#taken-in-m2) |
| D3 | Load/gen power on an absent phase is a model-level error | [Taken in M3, for loads and shunts](#taken-in-m3-for-loads-and-shunts) |
| D4 | Shunt at a node without all phases | [Taken in M3, for loads and shunts](#taken-in-m3-for-loads-and-shunts) |
| D5 | Absent-phase output is 0 everywhere | [Taken in M2](#taken-in-m2) |
| D6 | One exception class, `ReducedPhaseNotSupported` | [Taken in M2](#taken-in-m2) |
| D7 | Order of M3 and M4 | [Taken in M4](#taken-in-m4) |
| D8 | Identity rows in the Y-bus | [Taken in M3](#taken-in-m3) |
| D9 | Experimental gate | [Taken in M3](#taken-in-m3) |
| D10 | `PhaseConnection` numbering | [Taken in M3](#taken-in-m3) |
| D11 | Validation of `asym_line` input with `phases` | [Taken in M3](#taken-in-m3) |
| D12 | `WindingConnection1P` values | [Taken in M5](#taken-in-m5) |
| D13 | Polarity through `clock` | [Taken in M5](#taken-in-m5) |
| D14 | Earthing and the low-susceptance fallback | [Taken in M5](#taken-in-m5) |
| D15 | `u1`/`u2` are winding voltages | [Taken in M5](#taken-in-m5) |
| D16 | General phase shift | [Taken in M5](#taken-in-m5) |
| D17 | Loading from the complex sum per side | [Taken in M5](#taken-in-m5) |
| D18 | Working fixed taps | [Taken in M5](#taken-in-m5) |
| D19 | One centre-tapped value, to side only | [Taken in M6](#taken-in-m6) |
| D20 | Half-winding impedances `uk_half` and `pk_half` | [Taken in M6](#taken-in-m6) |
| D21 | `u2`, `uk`, `pk` and `sn` describe the full winding | [Taken in M6](#taken-in-m6) |
| D22 | Earthing of the centre tap | [Taken in M6](#taken-in-m6) |
| D23 | Loading of a centre-tapped side | [Taken in M6](#taken-in-m6) |
| D24 | Taps, clock, status, magnetising and faults of the centre-tapped model | [Taken in M6](#taken-in-m6) |
| D25 | Slot patterns for the start values | [Taken in M6, for start values](#taken-in-m6-for-start-values) |
| D26 | `LoadGenConnection` enum and default | [Taken in M7](#taken-in-m7) |
| D27 | Meaning of the delta entries | [Taken in M7](#taken-in-m7) |
| D28 | Absent phases for delta loads | [Taken in M7](#taken-in-m7) |
| D29 | Voltage dependence and per unit of delta loads | [Taken in M7](#taken-in-m7) |
| D30 | Per-phase output of delta loads | [Taken in M7](#taken-in-m7) |
| D31 | Short circuit with delta loads | [Taken in M7](#taken-in-m7) |
| D32 | Where delta is allowed | [Taken in M7](#taken-in-m7) |
| D33 | Calculation gates for delta loads | [Taken in M7](#taken-in-m7) |
| D34 | Model-level check of the fault phases | [Taken in M8](#taken-in-m8) |
| D35 | `InvalidShortCircuitPhases` with a new constructor | [Taken in M8](#taken-in-m8) |
| D36 | Default fault phases are not adapted to the node | [Taken in M8](#taken-in-m8) |
| D37 | The batch cache run swallows `InvalidShortCircuitPhases` | [Taken in M8](#taken-in-m8) |
| D38 | What Python checks and what stays in C++ | [Taken in M9](#taken-in-m9) |
| D39 | `asym_line.phases` validation | [Taken in M9](#taken-in-m9) |
| D40 | Sign of the mutual capacitance | [Taken in M9](#taken-in-m9) |
| D41 | `validate_single_phase_transformer` | [Taken in M9](#taken-in-m9) |
| D42 | `ValidationWarning` for the `u1`/`u2` rule | [Taken in M9](#taken-in-m9) |
| D43 | Load `connection` | [Taken in M9](#taken-in-m9) |
| D44 | Sensors on a `single_phase_transformer` | [Taken in M9](#taken-in-m9) |
| D45 | The experimental gate stays | [Taken in M9](#taken-in-m9) |
| D46 | User-manual page, no notebook | [Taken in M9](#taken-in-m9) |
| D47 | Maximum demand 20 %, minimum 10 %, power factor 0.95 | [Taken in M10](#taken-in-m10) |
| D48 | ISO coincidence factor at maximum demand only | [Taken in M10](#taken-in-m10) |
| D49 | Voltage regulators are fixed-neutral `YNyn0` banks | [Taken in M10](#taken-in-m10) |
| D50 | ISOs are two-winding | [Taken in M10](#taken-in-m10) |
| D51 | Split-phase transformers with an explicit Kersting split | [Taken in M10](#taken-in-m10) |
| D52 | The 44 phase-to-earth units of KLO14 are 12.7 kV / 250 V | [Taken in M10](#taken-in-m10) |
| D53 | $i_0 = p_0 = 0$ everywhere | [Taken in M10](#taken-in-m10) |
| D54 | Near-zero connectors are merged, zone-substation $u_k$ 1e-4 | [Taken in M10](#taken-in-m10) |
| D55 | References without `line`, short-circuit `node` and `asym_line` output | [Taken in M10](#taken-in-m10) |
| D56 | Tolerances of the real-network cases | [Taken in M10](#taken-in-m10) |
| D57 | Licence sidecars `MPL-2.0 AND BSD-3-Clause` | [Taken in M10](#taken-in-m10) |
| D58 | OpenDSS reference: $Z_1$/$C_1$ of one-phase lines, tolerance 1e-11 | [Taken in M10](#taken-in-m10) |
| D59 | Runtime cost in the debug test run | [Taken in M10](#taken-in-m10) |
| D60 | Per-method solver options in the validation parameters | [Taken in M11](#taken-in-m11) |
| D61 | Iterative current changes gated on a bus with an absent phase | [Taken in M11](#taken-in-m11) |
| D62 | KLO14 iterative current: IC-only `atol` of 1e-2 on the branch powers (to be reviewed by Michael) | [Taken in M11](#taken-in-m11) |
| D63 | Iterative current cases run with error tolerance 1e-10 and tighter tolerances | [Taken in M11](#taken-in-m11) |
| D64 | A NaN deviation raises `IterationDiverge` | [Taken in M11](#taken-in-m11) |
| D65 | Iterative current ignores `calculation_initialization` | [Taken in M11](#taken-in-m11) |
| D66 | Scope of iterative current in the validation cases | [Taken in M11](#taken-in-m11) |

D5 is also described in [Absent-phase output](#absent-phase-output), and D9 and D45 in the user manual page
[Reduced-phase networks](../user_manual/reduced-phase-networks.md#experimental-feature).

## Component models (phase domain)

### Per-unit convention for reduced-phase nodes

Keep `u_rated` as $\sqrt{3}$ × nominal phase-to-earth voltage for **every** node. The asymmetric per-unit base is
already phase-to-earth (`base_power_1p = 1e6/3` and `u_scale<asym> = 1/√3`, `<core>/common/common.hpp:93-96`).

| System | Phase-to-earth | `u_rated` |
|--------|---------------:|----------:|
| 12.7 kV SWER | 12.7 kV | 22 kV |
| 19.1 kV SWER | 19.1 kV | 33 kV |
| 240 V single-phase | 240 V | 415 V |
| 240-0-240 split-phase (leg to centre tap) | 240 V | 415.7 V |

This needs no new attribute, `ConflictVoltage` checks keep working, and it matches the Australian naming
convention. It must be stated clearly in the Node section of `docs/user_manual/components.md`.

### `asym_line` with `phases`

* **New attribute.** `phases`, of a new enum `PhaseConnection {abc = 0, a = 1, b = 2, c = 3, ab = 4, ac = 5,
  bc = 6}` (`IntS`, input-only, default `abc`; decision D10). Absent-phase impedances are left as NaN (decision
  D11).
* **Series admittance.** Let $P$ be the set of present phases. Fill the absent rows and columns of $Z$ with the
  identity, invert, and then zero the absent rows and columns:

  $$
  Y_{abc} = M_P \, \mathrm{inv}\!\left(Z_{PP} \oplus I\right) M_P
  $$

  Here $M_P$ is the diagonal 0/1 mask. The result is exact, because the absent block is decoupled. With a neutral,
  the absent phase conductors of the 4×4 matrix are replaced by an identity block before `kron_reduction`
  (`mask_absent_phase_conductors`). The helper `masked_inv(tensor, mask)` lives in
  `<core>/component/line_utils.hpp`. The neutral of a one- or two-phase line is present when `r_nn` or `x_nn` is
  given (for `abc` the existing rule, based on `r_na`/`x_na`, is unchanged).
* **Full-phase path.** For `abc` the constructor keeps the literal existing expressions, so three-phase results stay
  bit for bit identical.
* **Shunt.** Treat the capacitance matrix the same way. The `c0`/`c1` input form is only valid for `abc`.
* **Single end connected** (`asym_line.hpp:105-116`). The open-end branch shunt
  $\tfrac12 Y_{sh} + (Z + 2 Y_{sh}^{-1})^{-1}$ uses `masked_inv`, and the
  `(cabs(y_shunt) >= tol).all()` guard is applied to the present entries only. For `abc` this guard drops the whole
  open-end shunt when a *mutual* capacitance is 0; that pre-existing behaviour is kept and is a candidate for a
  separate fix, because changing it changes three-phase results. Michael decided on 2026-10-15 to leave it as it is
  for now.
* **SWER line.** A SWER line is `phases = a` (or `b` or `c`) with a 1×1 impedance whose self term already includes
  the Carson earth-return correction:

  $$
  r_{aa} = R_c + \pi^2 f \cdot 10^{-4}\ \Omega/\mathrm{km}, \qquad
  x_{aa} = 4\pi f \cdot 10^{-4} \ln\frac{D_e}{GMR}\ \Omega/\mathrm{km}, \qquad
  D_e \approx 658.5\sqrt{\rho / f}\ \mathrm{m}
  $$

  No neutral or Kron reduction is needed. This matches an OpenDSS `Line` with `phases=1` and `bus1=X.1`.
* **Symmetric parameters.** `sym_calc_param()` throws `ReducedPhaseNotSupported` for `phases != abc` (see
  [Calculation behaviour](#calculation-behaviour)).
* **Phase incidence.** `AsymLine::phase_incidence()` reports `phases → phases`, per phase; for `abc` this equals the
  default, so a model without reduced-phase lines keeps the topology fast path.

### New component: `single_phase_transformer`

A new branch component is better than extending `transformer`. The existing model is sequence-based, defined by
clock number, and enforces three-phase winding types.

**Input** (base `BranchInput`):

* `u1`, `u2`: the rated *winding* voltages, across the connection terminals (decision D15). `sn`, `uk`, `pk`,
  `i0`, `p0` as for `transformer`. There are no zero-sequence or winding-type fields; they have no meaning for a
  single-phase transformer.
* Tap fields, mirroring `transformer`, all optional (decision D18).
* `from_connection` and `to_connection`, of a new enum `WindingConnection1P`
  `{a = 1, b = 2, c = 3, ab = 4, ac = 5, bc = 6}` (decision D12), plus `ab_center_tapped = 7` for a centre-tapped
  (split-phase) winding, on the to side only (M6, decision D19).
* `uk_half`, `pk_half`: the short-circuit parameters between the primary and one half of a centre-tapped winding,
  optional, only for a centre-tapped to side (decision D20).
* `clock`: 0 or 6 (modulo 12); 6 reverses the polarity of the to winding (decision D13).
* `r_grounding_from`, `x_grounding_from`, `r_grounding_to`, `x_grounding_to`: the impedance between the earthed
  winding terminal and true earth, only for a phase-to-earth winding (decision D14). For SWER this is the earth
  electrode, which carries the entire SWER load current.

**Two-winding model.** Build the primitive admittance of the two coupled windings and map it onto the bus slots
with an incidence matrix $C$ (2×6). Row 1 of $C$ is $e_a - e_b$ for an `ab` winding, $e_a$ for an `a` winding, and
so on.

$$
Y_{prim} = y_s \begin{bmatrix} 1/n^2 & -1/n \\ -1/n & 1 \end{bmatrix} + Y_{mag}, \qquad
Y_{node} = C^{T} Y_{prim} C \in \mathbb{C}^{6\times 6}
$$

* $Y_{prim}$ is in siemens: $y_s$ and $Y_{mag}$ are referred to side 2 as in `transformer`, with
  $n = u_1'/u_2'$ from the tap-adjusted winding voltages. The magnetising admittance is split half on each side
  ($Y_{mag} = \tfrac{y_m}{2}\,\mathrm{diag}(1/n^2, 1)$), the same $\pi$ model as `transformer`. With one side
  open, $Y_{prim}$ is reduced to the open-end shunt $\tfrac{y_m}{2} + (1/y_s + 2/y_m)^{-1}$ of the connected
  winding (exactly 0 when $y_m$ is below `numerical_tolerance` in per unit).
* $Y_{node}$ is split into `yff`, `yft`, `ytf` and `ytt`, so `BranchCalcParam<asymmetric_t>` can hold it with no
  structural change. The slots outside the windings are set to exactly 0 by assignment at the end.
* Per-unit conversion: $Y_{pu}(s,t) = Y_{phys}(s,t)\, U_{b,s} U_{b,t} / S_{b,1p}$, where
  $U_b = u_{rated}/\sqrt3$ of each node.
* **Earthing impedance:** add an internal terminal for the earthed end of each phase-to-earth winding with a
  non-zero earthing impedance (an 8-terminal nodal matrix: 3 + 3 slots and 2 electrodes), stamp $1/z_g$ to earth,
  then eliminate the terminal with `kron_eliminate` (`<core>/component/line_utils.hpp`), a single-node Kron
  reduction on an Eigen matrix; eliminating the terminals one after the other equals eliminating them together.
* **Ungrounded windings:** a phase-to-phase winding with no earthed terminal has an undefined common-mode voltage.
  It gets the low-susceptance fallback of decision D14.

**Phase shift** (decision D16). `phase_shift()` ($\theta_{from} - \theta_{to}$) is

$$
\varphi(\text{to}) - \varphi(\text{from}) + \text{clock}\cdot\pi/6,
$$

where $\varphi$ is the angle of the winding voltage of a balanced positive-sequence set ($U_a = 1$,
$U_b = a^2$, $U_c = a$): $\varphi(a) = 0$, $\varphi(b) = -2\pi/3$, $\varphi(c) = 2\pi/3$, $\varphi(ab) = \pi/6$,
$\varphi(ac) = -\pi/6$ and $\varphi(bc) = -\pi/2$. For example:

* `ab` → `a`: $U_{ab} = \sqrt3 U_a \angle 30°$, so `phase_shift` = $-\pi/6$, and a SWER node behind it reports
  `u_angle` ≈ +30°;
* `bc` → `a`: $+\pi/2$;
* the design's `ca` → `a`, which is `ac` with clock 6: $\pi/6 + \pi \equiv -5\pi/6$;
* `a` → `b`: $-2\pi/3$, so the flat start puts slot b at $\theta_{from}$.

The earlier formulation (the winding angle relative to slot a, with the sign reversed) is the special case
to = `a`. The existing scalar phase-shift propagation in `tree_edge` (`<core>/topology.hpp:83`) and the M2 flat and
average-source starts then work unchanged, and output angles stay physical.

**Loading** (decision D17): $\max(|\sum_p s_{from,p}|, |\sum_p s_{to,p}|) / s_n$, with the complex sum of the
per-phase powers of each side.

**Calculations.** The incidence of a `single_phase_transformer` is never `abc` → `abc`, so the experimental gate
(D9), the symmetric routing and the state-estimation guard of M4 apply automatically; `sym_calc_param()` throws
`ReducedPhaseNotSupported` as a backstop. A winding side whose terminals are not all present at its node raises
`MissingPhaseForWinding` (D1).

**Validation against OpenDSS.** The OpenDSS references use a 1-phase, 2-winding wye-wye `Transformer` with
`kVs` = $u_1$, $u_2$, `ppm_antifloat=0` and the earth electrode as a `Reactor` (with $x > 0$) from node 4 of the bus
to earth; clock 6 swaps the two nodes of the to side. The cases use $i_0 = p_0 = 0$, because OpenDSS places the
magnetising branch differently from the $\pi$ model (the unit tests cover it), and from-side taps only. Lessons from
M5:

* The low susceptance of a phase-to-phase winding (decision D14) is visible in the results: without it, the
  transformer powers differ by about 5e-7 relative, and the short-circuit winding currents by up to about 1.5e-6 A.
  The reference therefore mirrors it as a 2-phase shunt `Reactor`, whose power and current are added to that side
  of the transformer. With it, PGM and OpenDSS agree to about 3e-10 relative, the same as for `asym_line`.
* OpenDSS cannot solve const_pq loads on a floating phase-to-phase secondary (an `ab` → `ab` winding, whose
  common mode is set by the loads only): its fixed-point iteration ends on the constant-impedance fallback of the
  loads, or diverges. The `phase-to-phase` case therefore uses mostly const_y loads there, with a small const_pq
  part.
* The iterative current method diverged on such a floating common mode, because its Y-bus holds the common mode only
  through the low susceptance and it treated the const_y loads as explicit currents. M5 therefore ran the
  `phase-to-phase` case with Newton-Raphson only, and `phase-to-phase-linear` (all const_y) with the linear methods.
  Since M11 the iterative current method puts the const_y loads into its matrix in a network with an absent phase
  (D61), which anchors the common mode, and `phase-to-phase` runs it too (11 iterations at error tolerance 1e-10).
  The limitation that remains is a floating common mode with only const_pq or const_i loads and no const_y load
  (see [Later work](#later-work)).
* OpenDSS replaces a fault resistance of 0 by 1e-4 Ω, so the reference models a bolted fault with 1e-13 Ω.
* The earth electrode carries the winding current of the isolating transformer, which the generator checks in
  OpenDSS and the short-circuit case compares with PGM. In a SWER fault it is the phasor sum of the fault current and
  the charging current of the SWER line, so its magnitude is close to, but not equal to, that of the fault current
  (188.816 A against 188.876 A for a bolted fault at the end of the SWER line).

**Centre-tapped (split-phase) secondary** (M6; decisions D19-D24). A `to_connection` of `ab_center_tapped` makes
the to winding a centre-tapped 240-0-240 V winding: half 1 from slot a to the centre tap $N$ and half 2 from $N$ to
slot b. The centre tap is earthed directly or through `r/x_grounding_to`. The legs map to **slots a and b**, both
measured to the earthed centre tap; slot c is absent. In normal operation $U_b \approx -U_a$; 240 V loads are wye
`asym_load` on a or b, and 480 V loads are delta a-b (M7).

* **Three-winding model.** The primary (0) and the halves (1, 2) are three coupled windings in the Dugan/OpenDSS
  N-winding form $Y = A^{T} Z_B^{-1} A$, with $A = \begin{bmatrix}1&-1&0\\1&0&-1\end{bmatrix}$ and
  $Z_B = z_0\mathbf{1} + z_b I$, in per unit on $s_n$ and the winding voltages $V_0 = u_1'$, $V_1 = V_2 = u_2'/2$
  (tap-adjusted as in the two-winding model).
  * $z_f = r + jx$ from `uk'`/`pk'` is the full winding (primary against both halves in series), and $z_{hs}$ the
    primary against one half, from `uk_half`/`pk_half` or the Kersting default $1.5r + j1.2x$ at the nominal tap
    (decision D20).
  * $z_b = 2(z_{hs} - z_{f,nom})$ is the branch of each half and $z_0 = z_f - z_b/2$ the branch of the primary, so
    $z_b + 2z_0 = 2z_f$ always.
  * In closed form: $Y_{00} = 1/z_f + y_m/2$, $Y_{01} = Y_{02} = -1/(2z_f)$,
    $Y_{11} = Y_{22} = (z_0 + z_b)/(2z_bz_f) + y_m/4$ and $Y_{12} = -z_0/(2z_bz_f)$. Without $y_m$ every row sums
    to 0. In siemens, $Y_{ij}\,s_n/(V_iV_j)$.
  * Two self-checks: with the primary shorted and the other half open, a leg sees $z_{hs}(u_2/2)^2/s_n$ (plus the
    earthing impedance of the centre tap); between a and b it sees $z_f u_2^2/s_n$. With $U_{to} = (v, -v, 0)$ the
    model gives the same currents as an `ab` secondary with the same `u2`.
* **Terminals.** The primary as in the two-winding model; half 1 is (+ slot a, − $N$) and half 2 (+ $N$, − slot b),
  where $N$ is earth when the centre tap is solidly earthed and the electrode terminal otherwise. Clock 6 swaps both
  pairs. The three windings are stamped on the 8-terminal nodal matrix, the electrode is eliminated with
  `kron_eliminate`, and the blocks are scaled to per unit as before; slot c is exactly 0 by assignment.
* **Phase shift.** $\varphi(\text{ab\_center\_tapped}) = 0$: the reference slot a is the positive terminal of half
  1. So `phase_shift` $= -\varphi(\text{from}) + \text{clock}\cdot\pi/6$: 0 for `a` → ct and $-\pi/6$ for `ab` → ct.
  The start values of slot b (180° from slot a, not 120°) come from the slot patterns of decision D25.
* **Status, magnetising, low susceptance, loading and taps** follow decisions D23 and D24: from side open gives the
  balancer $\frac{1}{2z_b}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$; the magnetising admittance is
  $\mathrm{diag}(y_m/2, y_m/4, y_m/4)$; a centre-tapped side has no low susceptance; each half is rated $s_n/2$;
  and the half branch is fixed at the nominal tap.

**Validation of the centre-tapped secondary against OpenDSS** (M6). The cases are generated by the same script as
for M5: power flow `tests/data/power_flow/single-phase-transformer/{split-phase,split-phase-flat-start,split-phase-linear}`
and short circuit `tests/data/short_circuit/single-phase-transformer-split-phase`. PGM and OpenDSS agree to about
3.4e-10 relative, as before; all cases pass with rtol 1e-8 and atol 1e-8 with a margin of at least 10x. `split-phase`
needed rtol 1e-6 in M6, because of the convergence criterion of iterative current at its default error tolerance of
1e-8; since M11 iterative current runs with an error tolerance of 1e-10 (D63) and `split-phase` passes at rtol 1e-8
(worst ratio 0.034). Lessons from M6:

* **Three-winding mapping.** A centre-tapped transformer is a 1-phase, 3-winding wye-wye-wye `Transformer` with
  `kVs` = $u_1$, $u_2/2$, $u_2/2$ and equal `kVAs` = $s_n$, so that all reactances are on $s_n$. The star
  resistances of the primary and half branches are given per winding (`%Rs` = $100\,\mathrm{Re}\,z_0$,
  $100\,\mathrm{Re}\,z_b$, $100\,\mathrm{Re}\,z_b$; `%loadloss` is not used with three windings), and the
  reactances between the windings as `XHL` = `XHT` = $100\,\mathrm{Im}\,z_{hs}$ and `XLT` = $200\,\mathrm{Im}\,z_b$.
  Half 1 is on `X.1.0` and half 2 on `X.0.2` (with an electrode, `X.1.4` and `X.4.2` plus a `Reactor` from `X.4` to
  earth); clock 6 reverses both pairs. The to-side values are summed per node over both halves.
* **YPrim check.** Before any network is solved, each centre-tapped transformer is checked on its own: the 6x6
  `CktElement.YPrim()` of OpenDSS must equal $C^{T} Y C$ of an independent numpy N-winding model
  ($A^{T} Z_B^{-1} A$, scaled to siemens) to 1e-9 relative (it agrees to 3e-16), and in open circuit
  $U_b = -U_a$ with $|U_a| = u_2/2$. This is the decisive check of the three-winding conventions (`%Rs` per winding,
  `XHL`/`XHT`/`XLT` on `kVAs[1]`, `kVs` and taps per winding), which were the highest risk of M6.
* **Centre-tap electrode.** The generator checks KCL at the centre tap in OpenDSS: the electrode current equals
  $I_a + I_b$, the phasor sum of the currents into the two halves at the legs. It carries the out-of-balance load
  current in power flow (19.02 A in `split-phase`) and the whole fault current of a leg-to-earth fault behind it
  (500.086 A); in a two-phase fault between the legs it carries no current.
* **`two_phase_to_ground`.** PGM joins the two phases and connects them to earth through $z_f$; bolted, both
  phases are at 0 V. The bolted fault is two OpenDSS `Fault` elements of 1e-13 Ω from each phase to earth, and the
  fault current per phase is the sum over all `Fault` elements of the fault and all their terminals. With
  $r_f > 0$, OpenDSS needs a bridge between the phases at a node with a non-zero voltage, and the large conductance
  of the bridge swamps the nodal admittance in double precision: for any bridge resistance the reference is only
  accurate to about 3e-9 relative at 240 V and 1e-8 at 22 kV, too coarse for these tolerances, so the generator
  does not produce it (PGM agrees with the best bridge values to about 1e-9). The same round-off makes a bolted
  two-phase fault (a 1e-13 Ω bridge) inexact at a low-impedance node in general: a reviewer measured 1.3e-4 relative
  at a three-phase 400 V node and 1.6e-2 at 22 kV. It is harmless only where the shorted phases sit near 0 V, as at a
  split-phase node, because then $\varepsilon\,|V|\,G$ (machine epsilon times the voltage times the 1e13 S bridge
  conductance) is negligible; the generator therefore solves every bolted two-phase fault again with a 10x larger
  bridge, which catches it, and requires the same fault current to 1e-10 relative. At a split-phase node a bolted `two_phase_to_ground` ab fault gives
  exactly the currents of the bolted `two_phase` fault: with both legs at 0 V the primary drives the halves antisymmetrically, so no
  current flows to earth at the centre tap. That scenario is therefore degenerate, and the case also has a bolted
  `two_phase_to_ground` bc fault at the 22 kV three-phase node 2, where the earth current is not zero and the
  fault differs from the two-phase fault (PGM and OpenDSS agree with a ratio of 0.035 at rtol/atol 1e-8).
* **Linear switch in the start tests.** A network with only constant-impedance loads, or none, switches power flow
  to the linear method (found in M6b; the rule is `all_const_y_` in `math_solver.hpp`). Cases and tests that exercise the slot patterns of
  the flat and average-source starts
  therefore need `const_pq` or `const_i` loads: `split-phase` and `split-phase-flat-start` have them, and
  `split-phase-linear` (all `const_y`) runs the linear methods.

### Phase-to-phase loads: `connection` on `asym_load` / `asym_gen`

Implemented in M7 (decisions D26 to D33).

* **Attribute.** `connection`, of the enum `LoadGenConnection {wye = 0, delta = 1}` (input-only, default `wye`;
  D26). For a delta connection the three `p_specified`/`q_specified` entries are the ab, bc and ca values (D27).
  A single phase-to-phase load is a delta load with one non-zero entry. The attribute is also present on
  `sym_load`/`sym_gen`, where `delta` is rejected.
* **Model.** With $D = \begin{bmatrix}1&-1&0\\0&1&-1\\-1&0&1\end{bmatrix}$, $U_\Delta = D U$ and, per active
  entry $k = (p, q)$ (ab = (a, b), bc = (b, c), ca = (c, a)):

  $$
  S_k(U) = S_k \left(\frac{|U_{\Delta,k}|}{U_{\Delta,nom}}\right)^n, \qquad
  I_{\Delta,k} = \overline{S_k(U) / U_{\Delta,k}}, \qquad I_{abc} = D^{T} I_\Delta, \qquad s = U \circ \overline{I}
  $$

  with $n = 0, 1, 2$ for `const_pq`, `const_i`, `const_y` and $U_{\Delta,nom}$ from the bus slot pattern
  (`slot_pattern(bus)`, D25): $\sqrt3$ for `three_phase` and $2$ for `split_phase` (D29). Equivalently, with
  $\alpha_p = U_p/U_{\Delta,k}$ and $\alpha_q = -U_q/U_{\Delta,k}$ ($\alpha_p + \alpha_q = 1$), entry k adds
  $g_r = S_k(U)\,\alpha_r$ to the phase powers of p and q, so $\sum_r g_r = S_k(U)$ (conservation).
* **Newton-Raphson.** The load power depends on the angles and magnitudes of the *other* phases at the same bus.
  With $\beta = S_k(U)\alpha_p\alpha_q$ and $\sigma = \begin{bmatrix}1&-1\\-1&1\end{bmatrix}$ over
  $r, s \in \{p, q\}$:

  $$
  \frac{\partial g_r}{\partial \theta_s} = j\sigma_{rs}\beta - n\,g_r\,\mathrm{Im}\,\alpha_s, \qquad
  V_s\frac{\partial g_r}{\partial V_s} = \sigma_{rs}\beta + n\,g_r\,\mathrm{Re}\,\alpha_s
  $$

  These go, with the solver convention $J = \partial PQ_{cal} - \partial PQ_{sp}$, into the H/N/M/L entries of the
  bus's own 3×3 diagonal block (`add_delta_load` in `newton_raphson_pf_solver.hpp`), so the sparsity pattern does not
  change. The rows satisfy two invariants that the tests use: $\sum_s \partial/\partial\theta_s = 0$ (rotation) and
  $\sum_s V_s \partial g_r/\partial V_s = n g_r$ (Euler homogeneity). The solver header documents the wye rule
  ($\partial PQ_{sp}/\partial\theta = 0$) and the delta terms.
* **Iterative current.** The right-hand side gets $D^T I_\Delta(U)$; the prefactorised matrix is unchanged.
* **Linear power flow and the Newton-Raphson linear start.** Constant admittance at nominal voltage:
  $Y_{add} = -D^T \mathrm{diag}(\overline{S_k}/U_{\Delta,nom}^2) D$, i.e. per active entry
  $y_k = -\overline{S_k}/U_{\Delta,nom}^2$ is added to $(p,p)$ and $(q,q)$ and subtracted from $(p,q)$ and $(q,p)$.
  Wye is the case $D = I$, $U_{nom} = 1$.
* **Results.** `calculate_load_gen_result` gives $i = D^T I_\Delta$ and $s = U \circ \overline{i}$ per phase (D30).
* **Exact zeros.** Inactive entries (exactly 0, or touching an absent phase) are skipped before any division, so
  $U_{\Delta,k} = 0$ on absent slots never produces NaN, and the absent slot rows and columns of the current,
  power, Jacobian and admittance are exactly +0. The existing absent-phase masks stay as the second layer. An active
  entry between two present phases at equal voltage would divide by 0; this is not physical in a power flow and
  surfaces as `IterationDiverge` or NaN.
* **Math (`<core>/math_solver/delta_load_gen.hpp`).** The pure functions `nominal_voltage`, `voltage_exponent`,
  `is_active`, `branch_power`, `phase_current`, `phase_power`, `power_and_jacobian` and `linear_admittance`,
  checked against closed forms (balanced three-phase and split-phase) and against central finite differences of
  `phase_power`.
* **Wiring.** `GenericLoadGen::connection()`; `ComponentTopology::load_gen_connection` (registered in
  `<core>/main_core/topology.hpp` only if a load_gen is delta, otherwise cleared), passed through
  `ReducedComponentTopology` and the supernode reduction, and coupled in `Topology::couple_all_appliance` into
  `MathModelTopology::load_gen_connection`, which a math model only gets (all `wye`, then `delta` by position) when
  one of its connected load_gens is delta (D32).
* **Alternative rejected.** A separate `delta_load` component would duplicate the load, sensor and output
  machinery.

**Validation of delta loads against OpenDSS** (M7). The cases are generated by the same script as for M5 and M6:
`tests/data/power_flow/asym-load-delta/{three-phase,three-phase-linear,two-phase,two-phase-linear,split-phase,split-phase-flat-start,split-phase-linear}`,
on an abc-only 11 kV network, on two-phase `ab` and `ac` laterals, and on the 240-0-240 V split-phase network of M6
with 480 V leg-to-leg loads. The references also hold the `asym_load` and `asym_gen` outputs. PGM and OpenDSS agree to
about 3.4e-10 relative, as before: `three-phase`, `three-phase-linear`, `split-phase-flat-start` and
`split-phase-linear` pass with rtol 1e-8 and atol 1e-8 with a margin of at least 10x. `two-phase` and `split-phase`
needed rtol 1e-6 in M7 (the convergence criterion of iterative current at its default error tolerance;
Newton-Raphson alone passed at 1e-8). Since M11 iterative current runs with an error tolerance of 1e-10 (D63):
`split-phase` passes at rtol 1e-8 and `two-phase` at rtol 1e-7 with atol 1e-8 (worst ratio 0.033).
`two-phase-linear` needs rtol 1e-7, and so does `two-phase` for the same reason, because the reactive power of the
lines nearly cancels there (6e-6 var on about 2 kvar, 2.8e-9 relative). Lessons from M7:

* **Delta loads in OpenDSS.** Each non-zero entry is a 1-phase `Load` with `conn=delta` and `bus1=X.p.q` (ca is
  `X.3.1`), with `kV` = the rated branch voltage: $u_{rated}$ at a node with the three-phase slot pattern, also at a
  two-phase `ab` or `ac` node, and $2u_{rated}/\sqrt3$ (480 V) at a split-phase node. Three equal entries can also be
  one native 3-phase `conn=delta` load with the total `kW`/`kvar`, which OpenDSS splits into a third per branch.
  `model` 1, 5 and 2 are const_pq, const_i and const_y, as for wye loads.
* **Base check.** After the solve, the generator checks every OpenDSS delta load: the sum of
  `CktElement.Powers()` over its conductors must equal $\sum_k S_k (|V_p - V_q| / kV)^n$ with the OpenDSS voltages
  to 1e-9 relative (it agrees to 6e-13). This is the decisive check of the conventions (the `kV` base of a 1-phase
  and a 3-phase delta load, and the voltage scaling of models 5 and 2), like the YPrim check of M6. It confirms that
  PGM's $U_{\Delta,nom}$ (D29) and OpenDSS's `kV` are the same base. It compares the element total, which for
  const_pq equals 3S whatever the split, so it does not prove the 3-phase split; that is proven end-to-end by the
  per-phase load outputs matching the OpenDSS terminal values.
* **Load outputs.** The per-phase `p`, `q` and `i` of a load are the sums per node of the terminal powers and
  currents of its OpenDSS `Load` elements (negated for a generator). They equal the PGM phase quantities of D30 to
  the same 3.4e-10, which confirms the phase-level output, including the reactive-power artefacts
  $\mp P/(2\sqrt3)$ of a resistive ab load at balanced voltages.
* **Linear switch.** As in M6, an all-`const_y` network switches to the linear method, so the Newton-Raphson and
  iterative current cases have const_pq and const_i delta loads; the `-linear` cases are all `const_y`.

**Validation of short-circuit phase validity against OpenDSS** (M8). The cases are generated by the same script as
for M5 to M7, in the `swer-chain` network without loads:
`tests/data/short_circuit/{swer-fault-impedance,swer-multiple-faults,swer-c-minimum}`, and the expected-error case
`short_circuit/swer-fault-on-absent-phase` (a fault on phase b at a SWER node raises `InvalidShortCircuitPhaseOrType`
with the fault ID, D34). PGM and OpenDSS agree to about 3.7e-10 relative, as before: all three pass with rtol 1e-8
and atol 1e-8 with a margin of at least 27x. All existing cases regenerate byte for byte. Lessons from M8:

* **Fault reactance in OpenDSS.** An OpenDSS `Fault` is resistive only, so a fault with $x_f \neq 0$ is a series
  R + jX `Reactor` (`bus2` = earth, or the second phase of a two-phase fault; $x_f > 0$ only, as for the
  electrodes). After the solve, the generator checks that the `Reactor` current equals $V / (r_f + jx_f)$ with the
  OpenDSS node voltages to 1e-9 relative (it agrees to 1e-16). This is the decisive check of the R/X convention,
  like the YPrim and base checks of M6 and M7.
* **Several faults per scenario.** A scenario is a tuple of simultaneous faults. The input has every fault ID, and
  the update switches the faults of a scenario on and the others off; a fault that is off has the null output
  ($i_f = 0$). All faults are single-phase to ground on a, as one calculation allows only one type and phase.
* **Upstream faults need an impedance.** A bolted fault at the start of a SWER line leaves the nodes behind it
  without voltage, so faults there get $i_f = 0$. That is correct, but it tests nothing, so in
  `swer-multiple-faults` every fault upstream of another one has an impedance, and every fault current is far
  from zero (129 A to 1.6 kA).
* **Voltage factor.** With minimum voltage scaling, *c* = 1.0 from the 22 kV source node, also for the fault at
  240 V; `swer-c-minimum` pins this (IEC 60909 would use 0.95 at the fault location). The reference scales the
  source voltage and impedance with the same *c* as PGM.

**Validation against real networks** (M10). The two networks `Network_1_Rural_SMR8` and `Network_2_Rural_KLO14` of
[Team-Nando/MV-LV-Networks](https://github.com/Team-Nando/MV-LV-Networks) (commit `8429f8a19f7b95acca2965cc48f2736880e5b31c`) have 3.4k
and 3.2k nodes, 10 and 7 SWER networks behind isolating transformers, voltage regulators and about 700 distribution transformers each. They test
the feature at a scale and with a data quality that the small synthetic cases do not have: a mix of
three-phase, split-phase and SWER transformers, and connectors with a near-zero impedance. The cases are generated by
`tests/data/generators/australian_swer_networks.py` (D47 to D59) and compared with OpenDSS: power flow in
`tests/data/power_flow/australian-swer/{smr8,klo14}` and short circuit in `tests/data/short_circuit/australian-swer/{smr8,klo14}`. PGM and OpenDSS
agree to about 1e-10 relative in node voltages and fault currents, and every case passes with a margin of at least 14x in the tolerance (worst
ratio 0.069 at M10; see D56 for the values after M11). The generated `README.md` of each case lists every deviation from the source repository's `backend.py`. The reduced-phase code needed no change. Lessons from M10:

* **Check that the load level is feasible before fixing it.** The 60 % of the transformer rating first approved does not converge
  in OpenDSS and collapses the voltage (D47). The planning step ran the networks in OpenDSS first, and the generator stops at a voltage gate of 0.8 pu.
* **OpenDSS one-phase lines use $Z_1$ and $C_1$ only**, and the matrix that a `LineCode` reports is not what the `Line` solves (D58).
  The OpenDSS tolerance had to be 1e-11, as 1e-12 stalls at round-off on SMR8.
* **The regularisation of PGM must be mirrored.** The zero-sequence low susceptance of a delta winding of a three-phase `transformer` is
  part of the PGM model, as is that of a phase-to-phase winding (D14). The reference therefore adds matching shunt `Reactor`s, as in M5 (D53).
* **Near-zero impedances limit the Newton-Raphson solver of PGM, with or without reduced-phase components** (D54). The three-phase-only version of each network stalls at
  the same deviation. An impedance floor, as an easy fix, changes the network (voltages by up to 0.4 %, losses by 4 % to 10 %); merging the connectors
  changes it by at most 5.5e-5 in voltage and 0.1 % in losses. Whichever fix is used, it must be applied to both tools,
  and the generator solves the unreduced network in OpenDSS to report how far the fix moved the solution. See [Later work](#later-work).
* **Size and runtime are design inputs.** The reference output was trimmed to what the budget allowed (D55: no three-phase `line`, short-circuit `node` or
  `asym_line` output; compact JSON; 12 significant digits), which keeps the data at 9.0 MB. The cases still add 136 s and 46 s to the debug test run (D59).
* **Short branches need a per-field `atol`** (D56): a power computed from the voltages at both ends of a branch of 1e-3 ohm
  carries the round-off of the voltages, so a purely relative tolerance on all fields would fail without a model error.

### Initial voltages

The `Vector(T)` rotation of 0°/−120°/+120° plus a scalar `phase_shift` is wrong for split-phase buses. It is fixed
with the `bus_slot_pattern` field (decision D25; the design called it `bus_slot_angle`):

* Each branch's incidence also reports the slot pattern of its to side (`to_pattern`): `split_phase` for a
  centre-tapped secondary, whose legs a and b are 180° apart, and `three_phase` for everything else. A SWER
  secondary needs no pattern: its only slot is a.
* The patterns spread from the centre-tapped windings along lines and links through buses without slot c, in the
  separate pass after the DFS, and are written together with the masks.
* The flat and average-source starts, and the iterative-current start, apply them only when `has_reduced_phase`
  is true.

### Sources

At first, sources are only allowed on nodes with all three phases. The balanced `y_ref` and the scalar `u_ref`
assume three phases. A source on a reduced-phase node raises an error.

## Calculation behaviour

Implemented in M4 (`<core>/main_model_impl.hpp`):

* **Flag.** `has_reduced_phase_components` is set at model construction when any branch reports a phase incidence
  other than abc → abc (`Branch::phase_incidence()`). That covers an `asym_line` with `phases != abc` (M3) and, once
  it exists, a `single_phase_transformer` (M5), without a per-type list. `phases` is input-only, so the flag never
  changes after construction.
* **Delta loads (M7).** A second flag, `has_delta_load_gens`, is set at construction when the component topology
  has a `load_gen_connection` vector, i.e. when any `asym_load` or `asym_gen` is delta (status and connectivity
  ignored; `connection` is input-only). The reduced-phase checks run first, so their message wins; then a model with
  delta load_gens raises `ReducedPhaseNotSupported` with the delta flavour for state estimation (both symmetries)
  and for symmetric power flow. The short-circuit symmetry selection is unchanged: loads are ignored in short
  circuits, so the symmetric short circuit stays valid (decision D31). The experimental gate checks the flag after
  the reduced-phase gate (decision D33).
* **Disconnected components count.** The flag looks at the input, not at the derived topology, so a reduced-phase
  line that is switched off (or a whole de-energised reduced-phase feeder) still marks the model. This is
  deterministic and independent of status updates in a batch.
* **State estimation** throws `ReducedPhaseNotSupported` (decision D6) for both symmetries and all methods.
* **Symmetric power flow** throws `ReducedPhaseNotSupported`.
* **Short circuit.** The symmetry selection is `is_three_phase && !has_reduced_phase_components`. A three-phase fault
  on a three-phase node of a mixed network runs through the asymmetric solver, which already supports
  `three_phase`. This also covers a short-circuit calculation without faults, which used to go to the symmetric
  solver.
* **Batch.** The check runs in `calculate()`, so it also runs in the batch cache run and is reported as one regular
  error, not as a batch error per scenario.
* **Backstops.** `MathSolver<symmetric_t>` (M2) and `AsymLine::sym_calc_param()` (M3) throw as well; they cannot be
  reached while the model-level check is in place. For delta load_gens (M7), `MathSolver` rejects state estimation
  and symmetric power flow when the math model has a `load_gen_connection` vector.

## Short-circuit specifics

* **Fault phase validity (M8, D34-D37).**
  * In `prepare_short_circuit_input` (`<core>/main_core/calculation_input_preparation.hpp`),
    `detail::check_fault_phases_present` checks, for every fault with status on in a math model with
    `has_reduced_phase`, that the resolved fault phases are a subset of the node's phase mask. Otherwise it throws
    `InvalidShortCircuitPhases` (`<core>/common/exception.hpp`) with the fault ID, the node ID and both phase sets,
    which maps to `InvalidShortCircuitPhaseOrType` in Python. Faults with status off and faults at isolated nodes are
    skipped (null output). The M2 check in the short-circuit solver, with the constructor
    `InvalidShortCircuitPhases(fault_phases, node_phases)`, stays as the backstop.
  * Default phases (single-phase faults on a, two-phase faults on bc; `Fault::get_fault_phase` in
    `<core>/component/fault.hpp`) are not adapted to the node (D36). They work on `a` SWER nodes and fail with a clear
    error, which says that the phases are the default, on `b`/`c` nodes and, for the two-phase faults, on
    split-phase and other nodes without phases b and c.
  * In a batch, the cache run swallows the error (D37), and every scenario checks its own faults with its own phase
    masks, so a scenario that fixes the fault succeeds and one that breaks it fails with the fault ID.
* **Valid cases:**
  * SWER phase-to-earth fault: `single_phase_to_ground` on a.
  * Split-phase leg-to-centre fault: `single_phase_to_ground` on a or b.
  * Split-phase leg-to-leg (480 V) fault: `two_phase` on ab.
  * Two-phase node: `single_phase_to_ground` on either phase, `two_phase` and `two_phase_to_ground` on both.
  * `three_phase` only at abc nodes.
* **One type per calculation.** The rule that all faults in one calculation share type and phase (`check_input_valid`)
  is more restrictive in mixed networks, e.g. a SWER on slot a combined with a fault on phase b elsewhere. Keep it,
  and document it (done in M8 in `docs/user_manual/components.md`).
* **Voltage factor *c*.** `Source::calc_param` picks *c* from the **source node** voltage level. IEC 60909 bases it
  on the nominal voltage at the fault location. This limitation exists today but matters more for long SWER and LV
  laterals. Documented in M8 (`docs/user_manual/components.md`); the M8 validation case `swer-c-minimum` pins the
  source-node *c*. Per-fault-location *c* is future work.
* **Earth return.** SWER phase-to-earth fault current is limited by:
  * the Carson self-impedance of the line, carried in `r_aa`/`x_aa`;
  * the isolating transformer's earth electrode (`r_grounding_to`);
  * any distribution-transformer earthing in the path;
  * `r_f`/`x_f`.

  Recommend minimum *c* with worst-case earthing resistances for protection reach studies. IEC 60909 does not
  cover earth-return systems specifically. Results follow the generic phase-domain model.

## Milestones

| # | Milestone | Main files | Tests | Effort / risk |
|---|-----------|-----------|-------|---------------|
| M0 | Fix `average_of_off_diagonal_of_matrix` | `<core>/common/matrix_utils.hpp` | Corrected expected values in `tests/cpp_unit_tests/component/test_asym_line.cpp` (which reproduced the bug); new `tests/cpp_unit_tests/common/test_matrix_utils.cpp` | XS / low |
| M1 | Phase-mask topology (mask always full; no behaviour change) | `<core>/calculation_parameters.hpp`, `<core>/common/exception.hpp`, `<core>/component/edge.hpp`, `<core>/component/transformer.hpp`, `<core>/component/generic_branch.hpp`, `<core>/main_core/topology.hpp`, `<core>/topology.hpp` | `tests/cpp_unit_tests/test_topology.cpp`: hand-built incidences, partial-winding error | M / medium |
| M2 | Solver masking (dangling-conductor and winding validation; Y-bus mask and identity; power-flow and short-circuit masking; state-estimation, symmetric and absent-phase-fault guards; output 0 on absent phases) | `<core>/topology.hpp`, `<core>/calculation_parameters.hpp`, `<core>/common/exception.hpp`, `<core>/component/link.hpp`, `<core>/math_solver/y_bus.hpp`, `newton_raphson_pf_solver.hpp`, `iterative_current_pf_solver.hpp`, `common_solver_functions.hpp`, `short_circuit_solver.hpp`, `math_solver.hpp` | `tests/cpp_unit_tests/test_topology.cpp`; `math_solver/test_y_bus.cpp`; new `math_solver/test_math_solver_pf_reduced_phase.cpp` (hand-computed single- and two-phase cases, equivalence with a full-phase network, bit-identical check with an all-abc mask, guards); `math_solver/test_math_solver_sc.cpp`; `component/test_node.cpp`; all existing validation cases unchanged | M / **high** (core numerics) |
| M3 (implemented) | `asym_line.phases` | `code_generation/data/attribute_classes/input.json`, `<core>/common/enum.hpp`, `src/power_grid_model/_core/enum.py`, `<core>/component/asym_line.hpp`, `<core>/component/line_utils.hpp`; regenerate code | Unit tests per phase set; validation cases against OpenDSS (`test.dss`): `tests/data/power_flow/asym-line/{single-phase-a,two-phase-ab,two-phase-bc-neutral-opendss-kron,swer-a,reduced-phase-feeder,reduced-phase-feeder-linear}` and `tests/data/short_circuit/{asym-line-reduced-phase,asym-line-swer-a}`, all generated by `tests/data/generators/reduced_phase_opendss_cases.py` (compared with PGM by `reduced_phase_compare_pgm.py` in the same folder). In most cases the generator gives OpenDSS the neutral-reduced matrix computed with numpy (cross-checked against an OpenDSS `LineCode` with `Kron=yes`); in `two-phase-bc-neutral-opendss-kron` OpenDSS does the neutral Kron reduction itself from the full matrices including the neutral, so that reference is independent of PGM's closed-form Kron. Expected-error cases `power_flow/asym-line/reduced-phase-symmetric-not-supported` and `state_estimation/reduced-phase-not-supported` | M / medium |
| M4 (implemented) | Symmetric-calculation guards | `<core>/main_model_impl.hpp`, `<core>/common/exception.hpp`, `src/power_grid_model/_core/error_handling.py`, `errors.py` | Native API tests (`tests/native_api_tests/test_api_model.cpp`); Python error-mapping and end-to-end tests; the expected-error validation cases of M3 | S / low |
| M5 (implemented) | `single_phase_transformer` (two windings, earthing, phase shift) | New `<core>/component/single_phase_transformer.hpp`; `<core>/all_components.hpp`; `<core>/main_core/main_model_type.hpp`; `input.json` and `dataset_definitions.json` (input, update, output and short-circuit output as a branch); enum `WindingConnection1P` | `tests/cpp_unit_tests/component/test_single_phase_transformer.cpp` (compare with hand-computed $C^T Y C$ for every connection pair; earthing impedance, low susceptance, clock, status, taps, `kron_eliminate`); `tests/cpp_unit_tests/test_topology.cpp`; native API tests (`tests/native_api_tests/test_api_model.cpp`); Python error-mapping tests (`tests/unit/test_error_handling.py`); validation cases against OpenDSS, generated by `tests/data/generators/reduced_phase_opendss_cases.py`: power flow `tests/data/power_flow/single-phase-transformer/{swer-chain,phase-to-phase,phase-to-phase-linear}` (22 kV a-b → 12.7 kV SWER → 240 V, and phase-to-phase windings with clock 6), short circuit `tests/data/short_circuit/single-phase-transformer-swer` (SWER and 240 V faults, earth-electrode current, three-phase fault with M4 routing), and the expected-error case `power_flow/single-phase-transformer/missing-phase-for-winding` (see [the component](#new-component-single_phase_transformer)) | L / **high** |
| M6 (implemented) | Centre-tapped split-phase secondary | `<core>/component/single_phase_transformer.hpp` (`ab_center_tapped`, `uk_half`/`pk_half`), `<core>/common/enum.hpp`, `src/power_grid_model/_core/enum.py`, `input.json`; slot patterns (D25) in `<core>/calculation_parameters.hpp`, `<core>/topology.hpp`, `common_solver_functions.hpp`, `newton_raphson_pf_solver.hpp`, `iterative_current_pf_solver.hpp` | `tests/cpp_unit_tests/component/test_single_phase_transformer.cpp` (closed-form 3x3 primitive for every primary, Kersting defaults, overrides, equivalence with an `ab` secondary, electrode, clock 6, magnetising, status, taps, loading, invalid input); `tests/cpp_unit_tests/test_topology.cpp` and `math_solver/test_math_solver_pf_reduced_phase.cpp` (slot patterns); native API test with `max_iterations = 1` (`tests/native_api_tests/test_api_model.cpp`); Python tests (`tests/unit/test_error_handling.py`); validation cases against OpenDSS (`phases=1 windings=3`, checked by YPrim), generated by `tests/data/generators/reduced_phase_opendss_cases.py`: power flow `tests/data/power_flow/single-phase-transformer/{split-phase,split-phase-flat-start,split-phase-linear}` (unbalanced 240 V loads, solid and electrode-earthed centre taps, clock 6) and short circuit `tests/data/short_circuit/single-phase-transformer-split-phase` (leg-to-earth, leg-to-leg and two-phase-to-ground faults, fault return through the centre-tap electrode); 480 V loads follow with M7 (see [the component](#new-component-single_phase_transformer)) | M / medium-high |
| M7 (implemented) | Delta loads (`connection` on `asym_load`/`asym_gen`) | `input.json` (`connection`), `<core>/common/enum.hpp` and `src/power_grid_model/_core/enum.py` (`LoadGenConnection`), `<core>/component/load_gen.hpp`, `MathModelTopology`/`ComponentTopology` (`load_gen_connection`), new `<core>/math_solver/delta_load_gen.hpp`, `add_loads` in the Newton-Raphson and iterative current solvers, `common_solver_functions.hpp` (linear and results), `main_model_impl.hpp` (gates) | `tests/cpp_unit_tests/math_solver/test_delta_load_gen.cpp` (closed forms, finite-difference check of the Newton-Raphson Jacobian, exact zeros, linear admittance); `math_solver/test_math_solver_pf_reduced_phase.cpp` (KCL for all solvers, bit-identity with all load_gens flagged wye, quadratic convergence, backstops); `component/test_load_gen.cpp`, `test_topology.cpp`, `main_core/test_main_core_topology.cpp`, `common/test_exceptions.cpp`; native API tests (`tests/native_api_tests/test_api_model.cpp`); Python tests (`tests/unit/test_error_handling.py`); validation cases against OpenDSS (1-phase `Load conn=delta` on `bus.p.q` and a native 3-phase `conn=delta` load, checked by the branch-voltage base check), generated by `tests/data/generators/reduced_phase_opendss_cases.py`: `tests/data/power_flow/asym-load-delta/{three-phase,three-phase-linear,two-phase,two-phase-linear,split-phase,split-phase-flat-start,split-phase-linear}` (with `asym_load`/`asym_gen` outputs in the reference) and the expected-error cases `power_flow/asym-load-delta/{symmetric-not-supported,delta-on-absent-phase}` and `state_estimation/asym-load-delta-not-supported` (see [Phase-to-phase loads](#phase-to-phase-loads-connection-on-asym_load--asym_gen)) | M / medium |
| M8 (implemented) | Short-circuit phase validity | `<core>/main_core/calculation_input_preparation.hpp` (`check_fault_phases_present`), `<core>/component/fault.hpp` (`has_default_phase`), `<core>/common/exception.hpp` (`InvalidShortCircuitPhases` with the fault ID), `<core>/calculation_parameters.hpp` (`fault_phase_mask`), `<core>/job_adapter.hpp` (cache-run catch); decisions D34-D37 | `tests/cpp_unit_tests/common/test_exceptions.cpp`, `component/test_fault.cpp`, `math_solver/test_math_solver_sc.cpp` (solver backstop); native API tests (`tests/native_api_tests/test_api_model.cpp`: valid, absent, default and skipped faults); Python tests (`tests/unit/test_error_handling.py`: single and batch, per-scenario topology); validation cases against OpenDSS, generated by `tests/data/generators/reduced_phase_opendss_cases.py`: `tests/data/short_circuit/{swer-fault-impedance,swer-multiple-faults,swer-c-minimum}` (`x_f` as a `Reactor`, simultaneous faults on three-phase, SWER and 240 V nodes, minimum voltage scaling) and the expected-error case `short_circuit/swer-fault-on-absent-phase`; user docs in `docs/user_manual/components.md` | S / low |
| M9 (implemented) | Python validation and documentation (decisions D38-D46) | `src/power_grid_model/validation/_validation.py` (required values and ranges of `asym_line` with `phases`, new `validate_single_phase_transformer`, load `connection`: valid enum value, `delta` only on `asym_load` and `asym_gen`, no voltage regulator on a delta load_gen, sensor branch lists, the `ValidationWarning` hooks), `_rules.py`, `errors.py` (`NotApplicableValueError`, `SinglePhaseTransformerClockError`, `ValidationWarning`), `__init__.py`; the new user manual page `docs/user_manual/reduced-phase-networks.md` (a worked example instead of a notebook, D46); `docs/user_manual/{components,calculations,data-validator,model-validation}.md`, `docs/advanced_documentation/c-api.md` and `docs/algorithms/sc-algorithms.md`; the comment of the experimental gate in `<core>/main_model_impl.hpp` | `tests/unit/validation/test_reduced_phase_validation.py` (every rule, with the mirror tests against C++), `test_reduced_phase_cases.py` (all reduced-phase, split-phase and delta validation cases validate cleanly without a warning), `test_rules.py`, `test_errors.py`, `test_input_validation.py`, `test_validation_functions.py`; `tests/unit/test_docs_reduced_phase_example.py` (runs the code of the user manual page) | M / low |
| M10 (implemented) | Validation against two real Australian rural networks with SWER: `Network_1_Rural_SMR8` and `Network_2_Rural_KLO14` from [Team-Nando/MV-LV-Networks](https://github.com/Team-Nando/MV-LV-Networks) (commit `8429f8a19f7b95acca2965cc48f2736880e5b31c`, BSD-3-Clause, AusNet Services data). Only the 22 kV and SWER (12.7 kV) parts are modelled; the LV side is lumped into one constant-power load per distribution transformer (three-phase, centre-tapped 250-0-250 V, SWER 12.7 kV / 240 V, and the 44 phase-to-earth units of KLO14) on its secondary busbar. Minimum demand is 10 % and maximum demand 20 % of the transformer kVA, at power factor 0.95. The generator writes the OpenDSS script and the PGM input from the same parsed data, following the conventions of the source repository's `backend.py`; decisions D47-D59 | `tests/data/generators/australian_swer_networks.py` (new: downloads the workbooks at the pinned commit and checks their sha256, maps them to OpenDSS and PGM, runs the self-checks, runs OpenDSS and writes the cases); `tests/data/generators/reduced_phase_compare_pgm.py` (extended: compact JSON rows, power-flow batches, a wider tolerance grid when no standard candidate passes); `LICENSES/BSD-3-Clause.txt`; `tests/unit/validation/test_reduced_phase_cases.py` (`MIN_CASE_COUNT` 37); no library changes | Four case folders, each with a generated `README.md` that lists every deviation from `backend.py`: power flow `tests/data/power_flow/australian-swer/{smr8,klo14}` (asymmetric batch of two scenarios, minimum and maximum demand; `newton_raphson` on both, `iterative_current` on SMR8 only until M11) and short circuit `tests/data/short_circuit/australian-swer/{smr8,klo14}` (IEC 60909, maximum voltage scaling, a batch of 8 (SMR8) and 9 (KLO14) faults on 22 kV, SWER and LV busbars). Measured agreement with OpenDSS: power flow rtol 1e-6, atol 1e-8 (1e-3 on `[pq]_(from|to)`), worst ratio 0.066 (SMR8) and 0.069 (KLO14), node voltages to 1.1e-10 and 9.8e-11 relative; short circuit rtol 1e-8, atol 1e-8 (1e-6 on `i_(from|to)` and its angle), worst ratio 0.033 and 0.030, fault currents to 6.7e-11 and 9.2e-11 relative. The data is 4.9 MB for power flow and 4.1 MB for short circuit | M / medium (data size, convergence, conditioning) |
| M11 (implemented) | The iterative current method works for every reduced-phase, split-phase, delta and real-network power-flow case (decisions D60-D66). Solver: const_y loads in the iterative current matrix at networks with an absent phase, refactorised when they change; a NaN deviation raises `IterationDiverge`; NaN folding of the maximum deviation over buses and phases. Validation: `iterative_current` with `error_tolerance` 1e-10 on every case that lists `newton_raphson` (including `australian-swer/klo14`, `single-phase-transformer/phase-to-phase`, both `split-phase-flat-start` cases and the expected-error cases) and tighter tolerances. The residual (defect-correction) form C1 was tried and reverted (D62) | `<core>/math_solver/iterative_current_pf_solver.hpp` (`absent_phase_` gate, cached const_y admittance, refactorisation), `<core>/math_solver/iterative_pf_solver.hpp` (NaN check), `<core>/math_solver/common_solver_functions.hpp` (`add_linear_load`, `add_const_y_loads`, `fold_max_deviation`), `<core>/math_solver/newton_raphson_pf_solver.hpp`; `tests/cpp_validation_tests/test_validation.cpp` (per-method `extra_params`); `tests/data/generators/{reduced_phase_opendss_cases,australian_swer_networks,reduced_phase_compare_pgm}.py`; `params.json` and READMEs of the power-flow cases | `tests/cpp_unit_tests/math_solver/test_math_solver_pf_reduced_phase.cpp` (NaN deviation raises, NaN folding of any bus and phase, IC with const_y loads on a floating common mode, refactorisation when the const_y loads change, delta const_y loads at buses with an absent phase, linear current with exact const_y loads and a NaN load); validation cases with IC on all power-flow cases; a capture of 21506 result keys before and after each solver commit: NR, linear, symmetric and three-phase results bit-identical. Results: `phase-to-phase` IC 11 iterations (worst ratio 0.034, equal to Newton-Raphson); worst ratio against OpenDSS 0.033 to 0.080 on the M3 to M9 cases (D63), 0.044 on SMR8 and 0.018 on KLO14 for IC; the PF batch takes 256.6 s (was 238.6 s before KLO14 IC); ctest 323/323, pytest 2465 passed | M / medium-high (core numerics, gated) |
| Later | See [Later work](#later-work) | — | — | — |

### Cross-cutting concerns

* **Validation (Python):** implemented in M9 (decisions D38 to D46).
  * Python checks the data of each component and simple ID lookups:
    * `asym_line`: `phases` must be a valid enum value; every r, x and c entry between two present conductors is
      required; a value on an absent phase or a missing neutral, and `c0`/`c1` on a one- or two-phase line, are
      rejected (D39).
    * `single_phase_transformer`: required values and ranges, valid connections with centre-tapping on the to side
      only, `clock` 0 or 6, earthing only on a phase-to-earth winding, `uk_half`/`pk_half` only on a centre-tapped
      winding, and the tap limits (D41).
    * `connection` of the four load_gen types, `delta` only on `asym_load` and `asym_gen`, and no voltage regulator
      on a delta load_gen (D43).
    * Sensors on a `single_phase_transformer` are valid, a tap regulator on it is not (D44).
  * The warning-level rule of D15 and D21 is a `ValidationWarning` (D42).
  * Only C++ checks the rules that need the derived node phases: winding terminals, dangling conductors, load and
    shunt on absent phases, fault phases against node phases, and the calculation gates (D38). A valid result of the
    Python validation therefore does not guarantee that the calculation succeeds. A best-effort Python check is
    listed under [Later work](#later-work).
* **Code generation and C API:** automatic once the JSON definitions are updated (`uv run code_generation/code_gen.py`).
  New exceptions need a regex mapping in `src/power_grid_model/_core/error_handling.py`.
* **Batch updates:** all new attributes are input-only. Masks are rebuilt only when a topology-changing status
  update occurs.
* **Performance:** every new code path is gated on `has_reduced_phase` (per math model) or
  `has_reduced_phase_components` (per model). Benchmark asymmetric power flow on existing cases, for example
  `ieee14` and `matpower-case300`, to confirm there is no regression.
* **Main risks:**
  * Newton-Raphson convergence for delta and split-phase cases with a poor initial guess. Mitigated by the linear
    start and the slot patterns of decision D25.
  * Singularity from ungrounded secondaries. Mitigated by the low-susceptance fallback.
  * How faithfully the OpenDSS references model the centre-tap impedance split. Resolved in M6 by the YPrim
    check of every centre-tapped transformer.

## Later work

Everything below is outside the implemented milestones. Items that are also open questions were confirmed as out of
scope by Michael (2026-10-15).

* **State estimation with masks** (open question 5), with per-phase observability rules.
* **A tap regulator on a `single_phase_transformer`** (D18, D44).
* **Sources and voltage regulators on reduced-phase nodes** (open question 3).
* **A delta `sym_load` and symmetric delta calculations** (open question 2, D26, D33).
* **The voltage factor *c* per fault location** (open question 4).
* **An output attribute for the derived node phases**, to tell an absent phase from 0 V (see
  [Absent-phase output](#absent-phase-output)).
* **Branch-level delta output** (p_ab, p_bc, p_ca), which needs new solver output fields (D30).
* **A shared earth electrode** for the two windings of an `a` → `a` transformer (D14), and **loads to a floating
  neutral** (D22).
* **`ac` and `bc` centre-tapped windings** (D19) and **`uk_half_min/max`** (D24).
* **Start values behind a split-phase bus**, for a transformer fed from a split-phase node (D25).
* **The open-end shunt guard of `asym_line`**, which drops the whole open-end shunt when a mutual capacitance is 0
  (see [`asym_line` with `phases`](#asym_line-with-phases)); changing it changes three-phase results.
* **A best-effort Python check of the node phases** (D38).
* **A public Python experimental switch, or graduating the feature** (D45).
* **The M11 solver changes for all networks, upstream** (D61). The const_y loads in the IC matrix ("C2") are gated on a
  bus with an absent phase, so that three-phase results stay bit-identical. Applying them to every network makes IC converge
  in fewer iterations wherever const_y loads exist, but it changes three-phase results bit for bit, so it is a decision for the upstream maintainers.
* **The residual (defect-correction) form of IC, "C1"** (D62). Evaluated in M11 and not adopted: no gain in accuracy against
  OpenDSS or Newton-Raphson, and a stall floor of 2e-11 to 5e-11 in the deviation on SMR8 and KLO14 (see D62 for the numbers). It could be reconsidered
  only if a case shows a systematic IC offset that the form without C1 cannot remove. The patch is not in the repository.
* **An all-load preconditioner for IC at a floating common mode with only const_pq or const_i loads.** Since M11
  a floating common mode converges when some const_y load is at it. With only const_pq or const_i loads, IC stays
  finite but does not converge and raises `IterationDiverge` at `max_iterations`; the matrix could also hold the linearised
  admittance $-\overline{S}$ of those loads at the floating buses. OpenDSS cannot solve such a network either (M5).
* **The Newton-Raphson `SparseMatrixError` on KLO14 from the flat and average-source starts.** Both starts give a singular matrix on KLO14; the default (linear) start
  works, and IC (which always starts at the average source) converges. Not investigated.
* **IC honouring `calculation_initialization`** (D65). IC always starts from the average-source voltage.
* **The Newton-Raphson stall on networks with near-zero impedance branches** (D54). On the SMR8 and KLO14 data, with
  connectors down to 1.4e-7 ohm, the Newton-Raphson power flow stalls at a maximum deviation of about 3e-7, above the
  error tolerance of 1e-8. This is pre-existing: it also happens with only the three-phase parts. Possible remedies are
  merging zero-impedance branches automatically before the solve, or a convergence criterion relative to the conditioning of the admittance matrix. It
  changes the behaviour of three-phase networks too, so it needs a discussion with the upstream maintainers.

## Open questions

1. **Absent-phase output.** Resolved in M2: 0 everywhere, consistent with de-energised output. See
   [Absent-phase output](#absent-phase-output).
2. **Delta loads in symmetric calculations.** *Out of scope for now (confirmed by Michael, 2026-10-15); stays in the
   "Later" milestone row.* A delta load has a well-defined positive-sequence equivalent. Should
   symmetric calculations allow it?
3. **Sources on SWER or split-phase nodes.** *Out of scope for now (confirmed by Michael, 2026-10-15); stays in the
   "Later" milestone row.* Needed, for example, for islanded single-phase studies?
4. **Voltage factor *c*.** *Out of scope for now (confirmed by Michael, 2026-10-15); stays in the "Later" milestone
   row.* Is *c* per fault location required for protection studies on SWER and LV laterals?
5. **State estimation.** *Out of scope for now (confirmed by Michael, 2026-10-15); stays in the "Later" milestone
   row.* Is support with phase masks required? If so, what are the per-phase observability rules?
