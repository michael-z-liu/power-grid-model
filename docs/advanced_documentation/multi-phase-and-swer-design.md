<!--
SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>

SPDX-License-Identifier: MPL-2.0
-->

# Design proposal: reduced-phase components and SWER networks

```{note}
This is a design and investigation document. It describes what is needed to support one- and two-phase
components in power flow and short-circuit calculations. None of it is implemented yet.
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
  * The open-end shunt is only applied when *every* element of $Y_{shunt}$ is non-zero (line 110).
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

* **New fields.** Add to `MathModelTopology` (`<core>/calculation_parameters.hpp:165`):
  * `std::vector<std::uint8_t> bus_phase_mask`: bit 0 for a, bit 1 for b, bit 2 for c.
  * `bool has_reduced_phase`: a fast-path flag. When it is false the mask stays empty and nothing below runs.
  * `std::vector<std::array<double, 3>> bus_slot_angle`: see [Initial voltages](#initial-voltages).
* **Incidence per branch.** Add a virtual `Edge::phase_incidence()` that returns which slots are connected on the
  from and to side, and whether the branch is a winding. It defaults to abc → abc for all existing branches.
  Register it in `ComponentConnections` (`<core>/main_core/topology.hpp:261-300`).
* **Propagation.** Add `propagate_phases()` to `Topology::build_topology()` (`<core>/topology.hpp:139`), after
  `dfs_search()`. Seed each math model's source buses with `0b111`, then propagate breadth-first:
  * **Line-type branches:** a present slot on one side makes the same slot present on the other side.
  * **Transformer windings:** the far side becomes present only if *all* non-earth terminals of the feeding
    winding are present. Otherwise a winding could back-feed an absent slot and create a floating, singular
    sub-problem. A winding with only some terminals present raises a new `MissingPhaseForWinding` error.
  * **Links and supernodes:** take the union of the masks.
* **Rebuilds.** `phases` and the connection attributes are input-only. The mask therefore changes only through
  status updates, which already trigger a topology rebuild.

### Solver changes

| Area | Change |
|------|--------|
| Y-bus (`<core>/math_solver/y_bus.hpp:444`, `update_admittance_entries`) | After a diagonal block is summed, add `1.0` at `(p, p)` for each absent `p`. Branch parameters are exactly zero on absent slots, so that row and column are otherwise empty. The incremental update path reuses the same function. |
| Newton-Raphson (`newton_raphson_pf_solver.hpp`) | New `apply_absent_phase_constraints()`, modelled on `apply_pv_constraints()` (line 585). Set `H(p,p) = L(p,p) = 1`, zero the rest of the absent rows, and set the mismatch to 0. The update keeps $V_p = 0$ and never divides by it. |
| Newton-Raphson starts (`make_flat_start`, `make_linear_start`, lines 289-339) | The flat start multiplies by the mask and uses `bus_slot_angle`. The linear start works as is, because of the identity in Y. |
| Iterative current (`iterative_current_pf_solver.hpp:172-195`), `iterative_pf_solver.hpp:110-117` | Mask the load right-hand side, which would otherwise be `0/0`, and the initial voltage. |
| Linear power flow | No change; the identity in Y gives $U_p = 0$. |
| Results (`common_solver_functions.hpp`) | Mask `conj(s/u)` to 0 on absent slots. |
| Short circuit | The Y-bus identity is enough. Fault substitution only touches fault phases, which are validated to be present (see [Short-circuit specifics](#short-circuit-specifics)). |
| State estimation | Not supported at first: throw `NotImplementedError` when `has_reduced_phase`. The observability check works per bus, not per phase, and pivot perturbation would hide modelling errors. |
| Voltage regulators / PV buses | Not supported on reduced-phase buses at first. |
| Outputs (`<core>/main_core/output.hpp`, `<core>/component/node.hpp`) | Absent-phase `u`, `u_pu` and `u_angle` become NaN. Currents and powers are naturally 0. See the [open questions](#open-questions). |

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

* **New attribute.** `phases`, of a new enum `PhaseConnection {abc = 0, ab, ac, bc, a, b, c}` (`IntS`,
  input-only, default `abc`). Absent-phase impedances are left as NaN.
* **Series admittance.** Let $P$ be the set of present phases. Fill the absent rows and columns of $Z$ with the
  identity, invert, and then zero the absent rows and columns:

  $$
  Y_{abc} = M_P \, \mathrm{inv}\!\left(Z_{PP} \oplus I\right) M_P
  $$

  Here $M_P$ is the diagonal 0/1 mask. The result is exact, because the absent block is decoupled. With a neutral,
  do the same on the 4×4 matrix before `kron_reduction`. A shared helper `masked_inv(tensor, mask)` would live in
  `<core>/component/line_utils.hpp`.
* **Shunt.** Treat the capacitance matrix the same way. The `c0`/`c1` input form is only valid for `abc`.
* **Single end connected** (`asym_line.hpp:105-116`). The open-end branch shunt
  $\tfrac12 Y_{sh} + (Z + 2 Y_{sh}^{-1})^{-1}$ must also use `masked_inv`, and the
  `(cabs(y_shunt) >= tol).all()` guard must become per-present-phase.
* **SWER line.** A SWER line is `phases = a` (or `b` or `c`) with a 1×1 impedance whose self term already includes
  the Carson earth-return correction:

  $$
  r_{aa} = R_c + \pi^2 f \cdot 10^{-4}\ \Omega/\mathrm{km}, \qquad
  x_{aa} = 4\pi f \cdot 10^{-4} \ln\frac{D_e}{GMR}\ \Omega/\mathrm{km}, \qquad
  D_e \approx 658.5\sqrt{\rho / f}\ \mathrm{m}
  $$

  No neutral or Kron reduction is needed. This matches an OpenDSS `Line` with `phases=1` and `bus1=X.1`.
* **Symmetric parameters.** `sym_calc_param()` must throw for `phases != abc` (see
  [Calculation behaviour](#calculation-behaviour)).

### New component: `single_phase_transformer`

A new branch component is better than extending `transformer`. The existing model is sequence-based, defined by
clock number, and enforces three-phase winding types.

**Input** (base `BranchInput`):

* `u1`, `u2`, `sn`, `uk`, `pk`, `i0`, `p0`.
* Tap fields, mirroring `transformer`, so that a tap regulator can be added later.
* `from_connection` and `to_connection`, of a new enum `WindingConnection1P`:
  `{a, b, c, ab, bc, ca, ab_center_tapped, bc_center_tapped, ca_center_tapped}`. Centre-tapped values are allowed
  on the to side only.
* `r_grounding_from`, `x_grounding_from`, `r_grounding_to`, `x_grounding_to`: the impedance between the earthed
  winding terminal and true earth. For SWER this is the earth electrode, which carries the entire SWER load
  current.

**Two-winding model.** Build the primitive admittance of the two coupled windings and map it onto the bus slots
with an incidence matrix $C$ (2×6). Row 1 of $C$ is $e_a - e_b$ for an `ab` winding, $e_a$ for an `a` winding, and
so on.

$$
Y_{prim} = y_s \begin{bmatrix} 1/n^2 & -1/n \\ -1/n & 1 \end{bmatrix} + Y_{mag}, \qquad
Y_{node} = C^{T} Y_{prim} C \in \mathbb{C}^{6\times 6}
$$

* $Y_{node}$ is split into `yff`, `yft`, `ytf` and `ytt`, so `BranchCalcParam<asymmetric_t>` can hold it with no
  structural change.
* Per-unit conversion: $Y_{pu}(s,t) = Y_{phys}(s,t)\, U_{b,s} U_{b,t} / S_{b,1p}$, where
  $U_b = u_{rated}/\sqrt3$ of each node.
* **Earthing impedance:** add an internal neutral terminal, then Kron-reduce it away. This needs a generalised,
  dynamic-size `kron_reduction`.
* **Ungrounded windings:** a phase-to-phase winding with no earthed terminal has an undefined common-mode voltage.
  Reuse the low-susceptance fallback (`transformer.hpp:310-328`, `transformer_low_susceptance_ratio`).

**Phase shift.** `phase_shift()` ($\theta_{from} - \theta_{to}$) is the angle of the winding voltage relative to
slot a, with the sign reversed:

* `ab` winding: $U_{ab} = \sqrt3 U_a \angle 30°$, so `phase_shift` = $-\pi/6$;
* `bc` winding: $+\pi/2$;
* `ca` winding: $-5\pi/6$;
* reversed polarity adds $\pi$.

The existing scalar phase-shift propagation in `tree_edge` (`<core>/topology.hpp:83`) then stays correct for the
reference slot.

**Centre-tapped (split-phase) secondary.**

* Model it as a three-winding single-phase transformer: primary, half 1 and half 2, in the Dugan/OpenDSS
  N-winding form $Y_{prim} = A^{T} Z_B^{-1} A$.
* Half 1 is $U_a - U_N$ and half 2 is $U_N - U_b$. The centre tap $N$ is earthed directly or through the earthing
  impedance.
* Split-phase legs map to **slots a and b**, both measured to the earthed centre tap. Slot c is absent.
* In normal operation $U_b \approx -U_a$.
  * 240 V loads are wye `asym_load` on a or b.
  * 480 V loads are delta a-b.
* The half-winding impedances default to the usual approximations, with optional `uk_half`/`pk_half` overrides.
  Check them against Kersting, *Distribution System Modeling and Analysis*, and against OpenDSS.

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
* The DFS copies these patterns along with the phase shift.
* The flat and average-source starts use `bus_slot_angle` only when `has_reduced_phase` is true.

### Sources

At first, sources are only allowed on nodes with all three phases. The balanced `y_ref` and the scalar `u_ref`
assume three phases. A source on a reduced-phase node raises an error.

## Calculation behaviour

* **Symmetric power flow and state estimation** throw a new `SymmetricCalculationNotSupported` error (or
  `NotImplementedError`) if the model contains any of:
  * an `asym_line` with `phases != abc`;
  * a `single_phase_transformer`;
  * (optionally) a delta `asym_load` or `asym_gen`.

  `phases` is input-only, so the check is a flag cached at model construction. Add it next to the existing
  method checks in `<core>/main_model_impl.hpp`. `AsymLine::sym_calc_param()` also throws as a backstop.
* **Short circuit.** Change the symmetry selection at `main_model_impl.hpp:362-368` to
  `is_three_phase && !has_reduced_phase_components`. A three-phase fault on a three-phase node of a mixed network
  then runs through the asymmetric solver, which already supports `three_phase`.

## Short-circuit specifics

* **Fault phase validity.**
  * In `prepare_short_circuit_input` (`<core>/main_core/calculation_input_preparation.hpp`), check that the fault
    phases are a subset of the node's phase mask. Otherwise throw `InvalidShortCircuitPhases`
    (`<core>/common/exception.hpp:248`) with an extended message.
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
| M1 | Phase-mask topology (mask always full; no behaviour change) | `<core>/calculation_parameters.hpp`, `<core>/component/edge.hpp`, `<core>/main_core/topology.hpp`, `<core>/topology.hpp` | `tests/cpp_unit_tests/test_topology.cpp`: hand-built incidences, partial-winding error | M / medium |
| M2 | Solver masking | `<core>/math_solver/y_bus.hpp`, `newton_raphson_pf_solver.hpp`, `iterative_current_pf_solver.hpp`, `iterative_pf_solver.hpp`, `common_solver_functions.hpp`, `short_circuit_solver.hpp`, `math_solver.hpp` (state-estimation rejection), `<core>/main_core/output.hpp`, `<core>/component/node.hpp` | Math solver unit tests with a forced mask; all existing validation cases unchanged | M / **high** (core numerics) |
| M3 | `asym_line.phases` | `code_generation/data/attribute_classes/input.json`, `<core>/common/enum.hpp`, `src/power_grid_model/_core/enum.py`, `<core>/component/asym_line.hpp`, `<core>/component/line_utils.hpp`; regenerate code | Unit tests per phase set; validation cases `tests/data/power_flow/asym-line/{single-phase-a,two-phase-ab,swer-a}` with OpenDSS `test.dss`; short-circuit cases | M / medium |
| M4 | Symmetric-calculation guards | `<core>/main_model_impl.hpp`, `<core>/common/exception.hpp`, `src/power_grid_model/_core/error_handling.py`, `errors.py` | C++ main-model tests; Python error-mapping tests | S / low |
| M5 | `single_phase_transformer` (two windings, earthing, phase shift) | New `<core>/component/single_phase_transformer.hpp`; `<core>/all_components.hpp`; `<core>/main_core/main_model_type.hpp`; `input.json` and `dataset_definitions.json` (input, update, output and short-circuit output as a branch); enum `WindingConnection1P` | `test_single_phase_transformer.cpp` (compare with hand-computed $C^T Y C$); validation case 22 kV a-b → 12.7 kV SWER → 240 V against OpenDSS | L / **high** |
| M6 | Centre-tapped split-phase secondary | Same component; dynamic-size Kron reduction in `line_utils.hpp`; slot angles | OpenDSS centre-tap transformer (`phases=1 windings=3`, `bus2=X.1.0`, `bus3=X.0.2`) with unbalanced 240/480 V loads | M / medium-high |
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

1. **Absent-phase output.** Should absent-phase voltages be NaN (proposed) or 0? NaN is more honest; 0 matches the
   current convention for non-energised nodes.
2. **Delta loads in symmetric calculations.** A delta load has a well-defined positive-sequence equivalent. Should
   symmetric calculations allow it?
3. **Sources on SWER or split-phase nodes.** Needed, for example, for islanded single-phase studies?
4. **Voltage factor *c*.** Is *c* per fault location required for protection studies on SWER and LV laterals?
5. **State estimation.** Is support with phase masks required? If so, what are the per-phase observability rules?
