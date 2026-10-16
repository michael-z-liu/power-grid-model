<!--
SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>

SPDX-License-Identifier: MPL-2.0
-->

# Design proposal: reduced-phase components and SWER networks

```{note}
This is a design and investigation document. It describes what is needed to support one- and two-phase
components in power flow and short-circuit calculations. Milestones M0 to M5 are implemented. M3 makes one- and
two-phase lines reachable from the input data, as an experimental feature (decision D9).
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

### Pre-existing defect found during this investigation

`average_of_off_diagonal_of_matrix` (`<core>/common/matrix_utils.hpp:16-18`) adds `matrix(1, 2)` twice and leaves
out `matrix(0, 2)`. `AsymLine::sym_calc_param` uses it, so symmetric results for any non-transposed `asym_line`
are slightly wrong. This is fixed in milestone M0.

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
  * `bus_slot_angle` (see [Initial voltages](#initial-voltages)) is not part of the topology milestone and is added
    together with the start values (M6).
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
| Newton-Raphson starts | The flat, average-source and linear starts mask the initial voltage of absent slots to exactly 0. The linear start also masks the load admittance per bus. Slot angles (`bus_slot_angle`) come with the split-phase work (M6). |
| Iterative current (`iterative_current_pf_solver.hpp`) | Mask the initial voltage and, after the load and source loop, the right-hand side (which holds `conj(s/0)` = NaN on absent slots). The prefactorised matrix already has the identity rows. `linear_current` is covered as well. |
| Linear power flow (`common_solver_functions.hpp`, `add_linear_loads`) | The load admittance $-\overline{s}$ is masked per bus. Otherwise it would be added to the identity on absent slots and give $(1 - \overline{s_p}) U_p = 0$, which is singular for $s_p = 1$. |
| Results (`common_solver_functions.hpp`) | The load/gen `s` and `i` (`conj(s/u)`) are set to 0 on absent slots. Bus injection is exactly 0 there because $U_p = 0$. |
| Short circuit (`short_circuit_solver.hpp`) | The Y-bus identity carries through `copy_y_bus`. Fault substitution only touches fault phases. In a reduced-phase model the solver checks that the fault phases are present at the faulted bus and raises `InvalidShortCircuitPhases` otherwise, so a three-phase fault, or a fault on an absent phase, at a reduced bus is rejected instead of silently giving zero fault current. The ID-based check in input preparation follows in M8. |
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

Taken in M2:

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

Taken in M4:

* **D7. Order of M3 and M4.** M4 lands before M3, so that reduced-phase models are never routed to the symmetric
  solver. A model with reduced-phase components always runs the asymmetric short circuit, also when every fault is
  `three_phase` and when there is no fault at all.

Taken in M3 (decided by recommendation; the points marked "decided by Michael" were reviewed on 2026-10-15):

* **D8. Identity rows in the Y-bus.** `YBus::admittance()` keeps the identity rows of absent slots; this is
  documented at the accessor. Consumers that need the physical admittance (future state estimation) must mask these
  slots.
* **D9. Experimental gate** (decided by Michael). Until the Python validation (M9)
  exists, a calculation on a model with one- or two-phase branches raises `ExperimentalFeature` ("A network with
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

Taken in M3, for loads and shunts:

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

Taken in M5, for `single_phase_transformer` (all decided by Michael, 2026-10-15):

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
  is "Later").
* **D15. `u1`/`u2` are winding voltages**, the rated voltages across the connection terminals: for example 22 kV
  for an `ab` primary, 12.7 kV for an `a` SWER secondary and 240 V for an `a` LV winding. Node `u_rated` stays
  $\sqrt3 \times$ the phase-to-earth voltage, so the nominal winding voltage implied by a node is $u_{rated}/\sqrt3$
  for a phase-to-earth and $u_{rated}$ for a phase-to-phase connection. Off-nominal ratios follow automatically and
  are legitimate, so C++ does not check `u1` against the node voltage; M9 adds a warning-level Python rule
  ($|u/U_{nom} - 1| > 0.2$) for the classic mistake of `u1 = 22 kV` on an `a` winding at a 22 kV node.
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
  otherwise the short circuit parameters stay at `uk`/`pk` on that side of `tap_nom`. The tap regulator stays
  "Later"; a `transformer_tap_regulator` on a `single_phase_transformer` raises `InvalidRegulatedObject`.

Taken in M6, for the centre-tapped (split-phase) secondary of `single_phase_transformer` (all decided by Michael,
2026-10-16):

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
  * **Valid ranges** (M9 Python rules; no C++ check, as `transformer` does not check `uk = 0` either):
    `uk_half >= pk_half / (sn/2)`; `uk_half > uk/2` (otherwise the half branch $z_b$ is 0 or negative, and singular
    at equality); `pk/4 <= pk_half <= pk/2` (non-negative branch resistances).
* **D21. `u2`, `uk`, `pk` and `sn` describe the full winding.** `u2` is the full winding, end to end: 480 V for
  240-0-240 V, consistent with D15 (the voltage across the connection terminals a-b). `uk`, `pk` and `sn` are the
  full-winding nameplate: the primary against both halves in series, at rated `sn`. Each half is $u_2/2$. The
  split-phase node has `u_rated` $= 240\sqrt3 = 415.692$ V, because each leg is 240 V to the earthed centre tap; the
  nominal winding voltage implied by the node is therefore $2u_{rated}/\sqrt3$. M9 adds this form to the warning
  rule of D15.
* **D22. Earthing of the centre tap.** Solid by default, or through `r_grounding_to`/`x_grounding_to`, the M5
  fields; not provided means solid. `r/x_grounding_from` follows the M5 rule for the primary. Caveat: PGM has no
  neutral node, so wye loads on the split-phase node are leg-to-*true-earth*. With an electrode, the whole
  out-of-balance load current therefore returns through the electrode. That is right for earth faults (the current
  to true earth returns through the electrode) and pessimistic for loads with a neutral conductor (MEN). Loads to a
  floating neutral are "Later".
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
* The iterative current method diverges on such a floating common mode with const_pq loads, because its Y-bus holds
  the common mode only through the low susceptance. The `phase-to-phase` case therefore runs Newton-Raphson only,
  and `phase-to-phase-linear` (all const_y) runs the linear methods.
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
  The start values of slot b (180° from slot a, not 120°) come with the slot patterns of M6b.
* **Status, magnetising, low susceptance, loading and taps** follow decisions D23 and D24: from side open gives the
  balancer $\frac{1}{2z_b}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$; the magnetising admittance is
  $\mathrm{diag}(y_m/2, y_m/4, y_m/4)$; a centre-tapped side has no low susceptance; each half is rated $s_n/2$;
  and the half branch is fixed at the nominal tap.

### Phase-to-phase loads: `connection` on `asym_load` / `asym_gen`

* **New attribute.** `connection`, of a new enum `LoadGenConnection {wye = 0, delta = 1}` (input-only, default
  `wye`). For a delta connection the three `p_specified`/`q_specified` entries are the ab, bc and ca values. A
  single phase-to-phase load is a delta load with one non-zero entry.
* **Model.** With $D = \begin{bmatrix}1&-1&0\\0&1&-1\\-1&0&1\end{bmatrix}$ and $U_\Delta = D U$:

  $$
  I_\Delta = \overline{S_\Delta / U_\Delta}, \qquad I_{abc} = D^{T} I_\Delta
  $$

* **Constant current and impedance.** These scale with $|U_\Delta| / U_{\Delta,nom}$. $U_{\Delta,nom}$ comes from
  the bus slot pattern: $\sqrt3$ for a 120° pattern and $2$ for a 180° split-phase pattern.
* **Newton-Raphson.** The load power now depends on the angles and magnitudes of the *other* phases at the same
  bus. The derivatives go into off-diagonal entries of the bus's own 3×3 H/N/M/L blocks, so the sparsity pattern
  does not change. The solver header's assumption that $\partial PQ_{sp} / \partial\theta = 0$ no longer holds and
  must be updated.
* **Wiring.** Add `load_gen_connection` to `MathModelTopology`, next to `load_gen_type`.
* **Alternative rejected.** A separate `delta_load` component would duplicate the load, sensor and output
  machinery.

### Initial voltages

The `Vector(T)` rotation of 0°/−120°/+120° plus a scalar `phase_shift` is wrong for split-phase buses. Fix it with
the `bus_slot_angle` field:

* Each branch's incidence also reports the slot pattern on the far side. For a centre-tapped secondary that is
  `{0, π, ·}`; for a SWER secondary it is `{0, ·, ·}`.
* The patterns are copied along with the masks, in the separate pass after the DFS.
* The flat and average-source starts use `bus_slot_angle` only when `has_reduced_phase` is true.

### Sources

At first, sources are only allowed on nodes with all three phases. The balanced `y_ref` and the scalar `u_ref`
assume three phases. A source on a reduced-phase node raises an error.

## Calculation behaviour

Implemented in M4 (`<core>/main_model_impl.hpp`):

* **Flag.** `has_reduced_phase_components` is set at model construction when any branch reports a phase incidence
  other than abc → abc (`Branch::phase_incidence()`). That covers an `asym_line` with `phases != abc` (M3) and, once
  it exists, a `single_phase_transformer` (M5), without a per-type list. Delta loads (M7) would need an extra term.
  `phases` is input-only, so the flag never changes after construction.
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
  reached while the model-level check is in place.

## Short-circuit specifics

* **Fault phase validity.**
  * In `prepare_short_circuit_input` (`<core>/main_core/calculation_input_preparation.hpp`), check that the fault
    phases are a subset of the node's phase mask. Otherwise throw `InvalidShortCircuitPhases`
    (`<core>/common/exception.hpp`) with an extended message. Since M2 the short-circuit solver already rejects
    fault phases that are absent at the faulted bus, with the constructor
    `InvalidShortCircuitPhases(fault_phases, node_phases)`; M8 adds the check with the fault ID.
  * Default phases (single-phase faults on a, two-phase faults on bc; `<core>/component/fault.hpp:171-194`)
    are unchanged. They work on SWER nodes and fail with a clear error on split-phase nodes.
* **Valid cases:**
  * SWER phase-to-earth fault: `single_phase_to_ground` on a.
  * Split-phase leg-to-centre fault: `single_phase_to_ground` on a or b.
  * Split-phase leg-to-leg (480 V) fault: `two_phase` on ab.
* **One type per calculation.** The rule that all faults in one calculation share type and phase (`check_input_valid`)
  is more restrictive in mixed networks, e.g. a SWER on slot a combined with a fault on phase b elsewhere. Keep it,
  and document it.
* **Voltage factor *c*.** `Source::calc_param` picks *c* from the **source node** voltage level. IEC 60909 bases it
  on the nominal voltage at the fault location. This limitation exists today but matters more for long SWER and LV
  laterals. Document it now; per-fault-location *c* is future work.
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
| M6 | Centre-tapped split-phase secondary | Same component; `kron_eliminate` in `line_utils.hpp` (added in M5); slot angles | OpenDSS centre-tap transformer (`phases=1 windings=3`, `bus2=X.1.0`, `bus3=X.0.2`) with unbalanced 240/480 V loads | M / medium-high |
| M7 | Delta loads | `input.json` (`connection`), `<core>/component/load_gen.hpp`, `MathModelTopology`, `add_loads` in the Newton-Raphson, iterative current and linear solvers, `common_solver_functions.hpp` | Finite-difference check of the Newton-Raphson Jacobian; OpenDSS `Load conn=delta` and two-phase `bus.1.2` loads | M / medium |
| M8 | Short-circuit phase validity | `<core>/main_core/calculation_input_preparation.hpp`, `<core>/component/fault.hpp` | `tests/data/short_circuit/swer-*` (phase-earth with `r_f` and earthing resistance) against OpenDSS fault study | S / low |
| M9 | Python validation and documentation | `src/power_grid_model/validation/_validation.py` (required values, `validate_asym_line`, new `validate_single_phase_transformer`, load `connection`, sensor branch lists), `docs/user_manual/components.md`, `docs/user_manual/calculations.md`, `docs/algorithms/sc-algorithms.md`, example notebook | `tests/unit/validation/*` | M / low |
| Later | State estimation with masks; tap regulator on `single_phase_transformer`; sources on reduced-phase nodes; delta `sym_load`; *c* per fault location | — | — | — |

### Cross-cutting concerns

* **Validation (Python):**
  * `phases` must be a valid enum value.
  * Present phases need their r/x/c entries; absent phases must be NaN.
  * `c0`/`c1` are only allowed with `abc`.
  * Transformer connections must be compatible, and centre-tapping is allowed on the to side only.
  * Fault phase versus node phases depends on derived topology. Only C++ can check it fully; Python can do a
    best-effort check.
* **Code generation and C API:** automatic once the JSON definitions are updated (`uv run code_generation/code_gen.py`).
  New exceptions need a regex mapping in `src/power_grid_model/_core/error_handling.py`.
* **Batch updates:** all new attributes are input-only. Masks are rebuilt only when a topology-changing status
  update occurs.
* **Performance:** every new code path is gated on `has_reduced_phase` (per math model) or
  `has_reduced_phase_components` (per model). Benchmark asymmetric power flow on existing cases, for example
  `ieee14` and `matpower-case300`, to confirm there is no regression.
* **Main risks:**
  * Newton-Raphson convergence for delta and split-phase cases with a poor initial guess. Mitigated by the linear
    start and slot angles.
  * Singularity from ungrounded secondaries. Mitigated by the low-susceptance fallback.
  * How faithfully the OpenDSS references model the centre-tap impedance split.

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
