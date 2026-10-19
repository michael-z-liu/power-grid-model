<!--
SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>

SPDX-License-Identifier: MPL-2.0
-->

# Reduced-phase networks (SWER and split-phase)

Many distribution networks are not fully three-phase.
Two-phase and single-phase laterals branch off three-phase feeders, single-wire earth return (SWER) systems carry a
single conductor with the earth as the return path, and customers are often supplied from a 240-0-240 V split-phase
(centre-tapped) secondary.
This page explains how to model such networks in power-grid-model, shows a complete example and lists what is and is
not supported.

## Experimental feature

```{warning}
One- and two-phase lines, single-phase transformers and delta-connected loads and generators are an experimental
feature.
A calculation on a model that contains any of them raises an
{py:class}`ExperimentalFeature <power_grid_model.errors.ExperimentalFeature>` error, unless the experimental features
are enabled for that calculation.
Experimental features have no stability guarantees: their behaviour and their interface may still change.
```

The experimental features are enabled per calculation:

* **C API:** set the option before calling `PGM_calculate`:

  ```c
  PGM_Options* opt = PGM_create_options(handle);
  PGM_set_symmetric(handle, opt, PGM_asymmetric);
  PGM_set_experimental_features(handle, opt, PGM_experimental_features_enabled);
  ```

* **Python:** there is no public switch yet.
  The public `calculate_power_flow` and `calculate_short_circuit` methods do not accept it, but the private methods
  `_calculate_power_flow` and `_calculate_short_circuit` of `PowerGridModel` do:

  ```python
  result = model._calculate_power_flow(symmetric=False, experimental_features="enabled")
  ```

  These methods are private and may change or disappear in any release.

Constructing the model is not gated: invalid input is reported by `PowerGridModel(input_data)` whether or not the
experimental features are enabled.

## Concepts

### Phases of a node

Nodes have no attribute for their phases.
The phases of a node are derived from the components connected to it:

* every node with a source has all three phases;
* an [`asym_line`](components.md#asym-line) with `phases` (for example `a` for a SWER line) carries only those phases
  to the node at its other end;
* a [`single_phase_transformer`](components.md#single-phase-transformer) gives the node behind it the phases of its
  to-side winding (`a` for a SWER secondary, `a` and `b` for a centre-tapped secondary), provided that all phases of its
  from-side winding are present.

The phases follow the status of the branches, so they are derived again after a status update.

### Rated voltage of a node

The `u_rated` of **every** node stays $\sqrt3$ times its nominal phase-to-earth voltage, also when the node has only
one or two phases:

| System                                    | Phase-to-earth | `u_rated` |
|-------------------------------------------|---------------:|----------:|
| 12.7 kV SWER                              |        12.7 kV |     22 kV |
| 19.1 kV SWER                              |        19.1 kV |     33 kV |
| 240 V single-phase                        |          240 V |     415 V |
| 240-0-240 V split-phase (leg to centre tap) |        240 V |   415.7 V |

The rated voltages of a `single_phase_transformer` (`u1`, `u2`) are, in contrast, the voltages across the terminals of
its windings: 22 kV for a winding between phases a and b of a 22 kV node, 12.7 kV for a SWER winding and 480 V for a
240-0-240 V centre-tapped winding.

### Output of absent phases

A phase that is not present at a node is *absent*.
Every output of an absent phase is exactly 0: the node `u`, `u_pu`, `u_angle`, `p` and `q`, and all branch, appliance
and fault outputs on that phase.
This is the same convention as for a de-energised component; an absent phase and a phase at 0 V can only be told apart
from the input.

## Building blocks

| What                                       | How                                                                                                                                         | Reference                                                                         |
|--------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------|
| one- and two-phase lines, SWER lines       | `asym_line` with `phases` set to `a`, `b`, `c`, `ab`, `ac` or `bc`, and only the r, x and c entries of the present conductors                | [Asym Line](components.md#asym-line)                                              |
| single-phase transformers                  | `single_phase_transformer` with `from_connection` and `to_connection` (`a`, `b`, `c`, `ab`, `ac`, `bc`); `clock` 0 or 6 for the polarity     | [Single-Phase Transformer](components.md#single-phase-transformer)                |
| earth electrodes                           | `r_grounding_<side>` and `x_grounding_<side>` of a phase-to-earth winding or of the centre tap                                              | [Single-Phase Transformer](components.md#single-phase-transformer)                |
| 240-0-240 V split-phase secondaries        | `to_connection` `ab_center_tapped`, optionally with `uk_half` and `pk_half`                                                                 | [Centre-tapped secondary](components.md#centre-tapped-split-phase-secondary)       |
| fixed taps of a single-phase transformer   | the optional tap attributes; each tap limit is optional                                                                                     | [Single-Phase Transformer](components.md#single-phase-transformer)                |
| automatic tap changing                     | a `transformer_tap_regulator` that regulates a `single_phase_transformer` (asymmetric power flow); the transformer needs `tap_min`, `tap_max` and `tap_size`  | [Single-phase transformer regulated by a tap regulator](components.md#single-phase-transformer-regulated-by-a-tap-regulator) |
| phase-to-earth loads and generators        | `asym_load` and `asym_gen` with power only on the phases present at the node                                                                | [Generic Load and Generator](components.md#generic-load-and-generator)             |
| phase-to-phase (for example 480 V) loads   | `asym_load` and `asym_gen` with `connection` `delta`; the entries are the ab, bc and ca branches                                            | [Generic Load and Generator](components.md#generic-load-and-generator)             |
| faults at a SWER or split-phase node       | `fault` with a `fault_phase` that is present at the node                                                                                    | [Fault](components.md#fault)                                                      |

## Worked example

The example below models a SWER system:

* node 1 is a 22 kV three-phase node with a source;
* an isolating transformer 20 (22 kV between phases a and b, 12.7 kV phase-to-earth, 200 kVA) feeds the SWER node 2;
  its earth electrode has a resistance of 5 Ω;
* a SWER line 30 on phase `a` runs from node 2 to node 3; its self impedance includes the earth return;
* a distribution transformer 21 (12.7 kV phase-to-earth, 240-0-240 V centre-tapped, 25 kVA) feeds the split-phase
  node 4;
* at node 4, load 40 draws 3 kW on leg a and 2 kW on leg b (240 V each), and load 41 draws 4 kW between the legs
  (480 V).

<!-- example-start -->

First, the input data.
All attributes of the reduced-phase components that are not provided stay empty (NaN), as `initialize_array` leaves
them.

```python
import numpy as np

from power_grid_model import ComponentType, DatasetType, LoadGenType, PowerGridModel, initialize_array
from power_grid_model.enum import LoadGenConnection, PhaseConnection, WindingConnection1P
from power_grid_model.errors import ExperimentalFeature
from power_grid_model.validation import validate_input_data

# nodes: 22 kV three-phase, 12.7 kV SWER (u_rated 22 kV) and 240-0-240 V split-phase (u_rated 240 * sqrt(3))
node = initialize_array(DatasetType.input, ComponentType.node, 4)
node["id"] = [1, 2, 3, 4]
node["u_rated"] = [22e3, 22e3, 22e3, 240 * np.sqrt(3)]

source = initialize_array(DatasetType.input, ComponentType.source, 1)
source["id"] = 10
source["node"] = 1
source["status"] = 1
source["u_ref"] = 1.0
source["sk"] = 1e8

# 20: isolating transformer, 22 kV between phases a and b -> 12.7 kV phase-to-earth (SWER)
# 21: distribution transformer, 12.7 kV phase-to-earth -> 240-0-240 V centre-tapped
transformer = initialize_array(DatasetType.input, ComponentType.single_phase_transformer, 2)
transformer["id"] = [20, 21]
transformer["from_node"] = [1, 3]
transformer["to_node"] = [2, 4]
transformer["from_status"] = 1
transformer["to_status"] = 1
transformer["from_connection"] = [WindingConnection1P.ab, WindingConnection1P.a]
transformer["to_connection"] = [WindingConnection1P.a, WindingConnection1P.ab_center_tapped]
transformer["u1"] = [22e3, 12.7e3]
transformer["u2"] = [12.7e3, 480.0]
transformer["sn"] = [200e3, 25e3]
transformer["uk"] = [0.05, 0.04]
transformer["pk"] = [2e3, 400.0]
transformer["i0"] = 0.01
transformer["p0"] = [300.0, 60.0]
transformer["r_grounding_to"] = [5.0, np.nan]  # earth electrode of the SWER system; the centre tap is solidly earthed

# SWER line: one conductor on phase a; r_aa and x_aa include the earth return
line = initialize_array(DatasetType.input, ComponentType.asym_line, 1)
line["id"] = 30
line["from_node"] = 2
line["to_node"] = 3
line["from_status"] = 1
line["to_status"] = 1
line["phases"] = PhaseConnection.a
line["r_aa"] = 6.0
line["x_aa"] = 8.0
line["c_aa"] = 5e-8

# 40: 240 V loads on the legs a and b; 41: a 480 V load between the legs
load = initialize_array(DatasetType.input, ComponentType.asym_load, 2)
load["id"] = [40, 41]
load["node"] = 4
load["status"] = 1
load["type"] = LoadGenType.const_power
load["connection"] = [LoadGenConnection.wye, LoadGenConnection.delta]
load["p_specified"] = [[3e3, 2e3, 0.0], [4e3, 0.0, 0.0]]  # wye: a, b, c; delta: ab, bc, ca
load["q_specified"] = 0.0

input_data = {
    ComponentType.node: node,
    ComponentType.source: source,
    ComponentType.single_phase_transformer: transformer,
    ComponentType.asym_line: line,
    ComponentType.asym_load: load,
}
```

Validate the data for an asymmetric calculation (see [Validation](#validation)).
The result is `None`: the data is valid.

```python
errors = validate_input_data(input_data, symmetric=False)
```

The public `calculate_power_flow` raises `ExperimentalFeature`, because the experimental features are not enabled.
The private `_calculate_power_flow` can enable them (see [Experimental feature](#experimental-feature)).

```python
model = PowerGridModel(input_data)

try:
    model.calculate_power_flow(symmetric=False)
except ExperimentalFeature as error:
    gated = error

result = model._calculate_power_flow(symmetric=False, experimental_features="enabled")

u = result[ComponentType.node]["u"]  # phase-to-earth voltages in V, one row per node
u_angle = np.degrees(result[ComponentType.node]["u_angle"])
```

The voltages and angles are about:

| node | `u` (V)                    | `u_angle` (°)           |
|------|----------------------------|-------------------------|
| 1    | 12702.6, 12700.5, 12701.7  | −0.01, −120.01, 120.00  |
| 2    | 12694.0, 0, 0              | 29.85, 0, 0             |
| 3    | 12690.2, 0, 0              | 29.82, 0, 0             |
| 4    | 238.28, 238.43, 0          | 29.05, −150.92, 0       |

* **Absent phases are 0.** The SWER nodes 2 and 3 only have phase a, and the split-phase node 4 only has the legs a and
  b, so all other entries are exactly 0.
* **The SWER nodes are at about +30°.** The isolating transformer is fed by the phase-to-phase voltage
  $U_a - U_b = \sqrt3 U_a \angle 30°$, so its secondary leads phase a of the 22 kV node by 30°.
* **Leg b is opposite to leg a.** The distribution transformer has no phase shift (`a` → `ab_center_tapped`), so leg a
  of node 4 follows the SWER voltage (+29°, a little less because of the voltage drop) and leg b is 180° from it, at
  about −151°.
* **The 480 V load appears on both legs.** The output of a delta load is per phase at its node.
  Between two legs that are 180° apart, the 4 kW of load 41 is split about evenly: `p` is about 2.0 kW on each leg and
  `q` is about ±0.6 var, because the legs are not exactly opposite.

### Automatic tap changing

A `transformer_tap_regulator` can regulate the tap position of a `single_phase_transformer`, with the tap changing
strategies of the [regulated power flow](calculations.md#power-flow-with-automatic-tap-changing).
In the example, the isolating transformer 20 gets a tap changer on its 22 kV side and a regulator that keeps the voltage
on its SWER side between 12.8 kV and 13.0 kV.
The transformer needs `tap_min`, `tap_max` and `tap_size` (see [Validation](#validation)).
`u_set` and `u_band` are in volts across the controlled winding, here 12.9 kV and 200 V, and the voltage that is
regulated is the voltage of the SWER node to earth.
The extra code continues the example above.

```python
from power_grid_model.enum import BranchSide, TapChangingStrategy

# 20: a tap changer on the 22 kV side, taps -4 to 4 of 1.25 % (275 V) each; a higher position lowers the voltage
regulated_transformer = transformer.copy()
regulated_transformer["tap_side"][0] = BranchSide.from_side
regulated_transformer["tap_nom"][0] = 0
regulated_transformer["tap_pos"][0] = 0
regulated_transformer["tap_min"][0] = -4
regulated_transformer["tap_max"][0] = 4
regulated_transformer["tap_size"][0] = 275.0

# 50: keeps the voltage of the SWER side of transformer 20 (to earth, in volts) between 12.8 kV and 13.0 kV
regulator = initialize_array(DatasetType.input, ComponentType.transformer_tap_regulator, 1)
regulator["id"] = 50
regulator["regulated_object"] = 20
regulator["status"] = 1
regulator["control_side"] = BranchSide.to_side
regulator["u_set"] = 12.9e3
regulator["u_band"] = 200.0

regulated_input = {
    **input_data,
    ComponentType.single_phase_transformer: regulated_transformer,
    ComponentType.transformer_tap_regulator: regulator,
}
regulated_errors = validate_input_data(regulated_input, symmetric=False)

regulated_result = PowerGridModel(regulated_input)._calculate_power_flow(
    symmetric=False,
    experimental_features="enabled",
    tap_changing_strategy=TapChangingStrategy.any_valid_tap,
)
tap_pos = regulated_result[ComponentType.transformer_tap_regulator]["tap_pos"]  # -1
u_swer = regulated_result[ComponentType.node]["u"][1, 0]  # V, about 12855
```

Without the regulator the voltage of node 2 is 12694 V.
The regulator lowers the tap position to −1, which raises the voltage of node 2 to about 12855 V, inside the band.
A symmetric calculation of this model raises `ReducedPhaseNotSupported`.

<!-- example-end -->

## Validation

Pass `symmetric=False` to `validate_input_data` and `validate_batch_data`: reduced-phase networks only support
asymmetric calculations, and the asymmetric validation also checks the zero-sequence and per-phase attributes of the
three-phase components.

### What the Python validation checks

The Python validation checks the attributes of each component and the references between components:

* **`asym_line`:** `phases` must be a valid {py:class}`PhaseConnection <power_grid_model.enum.PhaseConnection>`.
  Every r, x and c entry between two present conductors is required (`MissingValueError`), and an entry on an absent
  phase or a missing neutral must be left empty (`NotApplicableValueError`).
  `c0` and `c1` are only allowed with `abc`.
  The self values must be positive; the mutual capacitances may have any sign.
* **`single_phase_transformer`:** the required attributes and the ranges of `u1`, `u2`, `sn`, `uk`, `pk`, `i0`, `p0`
  and the tap-dependent short-circuit parameters, as for `transformer`;
  * both connections must be valid, and `ab_center_tapped` is only allowed on the to side (`InvalidValueError`);
  * `clock` must be 0 or 6 modulo 12 (`SinglePhaseTransformerClockError`);
  * `r_grounding_<side>` and `x_grounding_<side>` are not allowed on a phase-to-phase winding
    (`NotApplicableValueError`), and the earthing resistance must not be negative;
  * `uk_half` and `pk_half` are only allowed on a centre-tapped winding, with `uk_half > uk / 2`,
    `pk / 4 <= pk_half <= pk / 2` and `uk_half >= pk_half / (sn / 2)`, using the defaults for a value that is not
    provided;
  * `tap_pos` must lie within the tap limits that are provided, and `tap_nom` between them when both are provided;
  * power and current sensors can measure it;
* **`transformer_tap_regulator` on a `single_phase_transformer`:** `control_side` must be a
  {py:class}`BranchSide <power_grid_model.enum.BranchSide>` (`InvalidAssociatedEnumValueError`).
  The regulated transformer must provide `tap_min`, `tap_max` and `tap_size` (`MissingValueError`), whatever the status
  of the regulator: without a limit the range of the tap changer is not bounded, which the calculation rejects with an
  `AutomaticTapInputError`, and with a `tap_size` of 0 the regulator would change nothing.
* **`connection` of the load_gen components:** it must be a valid
  {py:class}`LoadGenConnection <power_grid_model.enum.LoadGenConnection>`, `delta` is only allowed for `asym_load` and
  `asym_gen` (`InvalidValueError`), and a `voltage_regulator` cannot regulate a delta load or generator
  (`InvalidVoltageRegulationError`).

The validation also emits a `ValidationWarning` when `u1` or `u2` of a `single_phase_transformer` differs by more than
20% from the nominal winding voltage implied by the `u_rated` of its node, which is $u_{\text{rated}} / \sqrt3$ for a
phase-to-earth winding, $u_{\text{rated}}$ for a phase-to-phase winding and $2 u_{\text{rated}} / \sqrt3$ for a
centre-tapped winding.
A typical mistake it catches is `u1` = 22 kV on an `a` winding at a 22 kV node, where 12.7 kV is meant.
The warning does not make the data invalid (see [Validation warnings](data-validator.md#validation-warnings)).

The validation also emits a `ValidationWarning` when the `u_set` of a `transformer_tap_regulator` of a
`single_phase_transformer` differs by more than 20% from the nominal voltage of the controlled winding, which is
$u_{\text{rated}} / \sqrt3$ of the node for a phase-to-earth winding and for one leg of a centre-tapped winding, and
$u_{\text{rated}}$ for a phase-to-phase winding.
`u_set` is in volts across the winding, not at the `u_rated` of the node: the warning catches a `u_set` of 22 kV for a
12.7 kV SWER winding, or 480 V for the legs of a split-phase winding.
For batch data every scenario is checked after its update, and the warning is emitted once per call.

The validation also emits a `ValidationWarning` when a phase-to-phase winding (`ab`, `ac` or `bc`) of a
`single_phase_transformer` has no earth reference.
Such a winding has an undefined common-mode voltage; only a very small fallback admittance holds it, so the result
depends on the loads alone (see [Supported calculations](#supported-calculations)).
The network of the winding, up to the next transformer (the nodes reached through lines, `asym_line`s and links), has an
earth reference if it holds a source, a phase-to-earth or centre-tapped `single_phase_transformer` winding, a `wye_n` or
`zigzag_n` winding of a `transformer` or `three_winding_transformer`, a `shunt` with non-zero `g0` or `b0`, or a wye
`const_impedance` load or generator with power.
Constant-power and constant-current loads, delta loads, the capacitance of lines, and loads in a short circuit do not
count.
Only windings that are connected, switched on and energised from a source are checked.
For batch data every scenario is checked after its update, and the warning is emitted once per call and names the
scenarios.
The check is best-effort: it does not see the phases of the nodes.

### What only the calculation checks

The rules that need the phases of the nodes are only checked by the calculation, because these phases are derived from
the topology and change with the status of the branches.
They raise an error from `power_grid_model.errors`; the errors about a single component name its ID:

| Rule                                                                                                             | Error                              |
|------------------------------------------------------------------------------------------------------------------|------------------------------------|
| all phases of a transformer winding must be present at its node                                                  | `MissingPhaseForWinding`           |
| a conductor of a line must be present at the nodes it connects                                                   | `DanglingConductor`                |
| a load or generator may not have power on an absent phase; a delta entry needs both of its phases                | `LoadGenOnAbsentPhase`             |
| a `sym_load` or `sym_gen` with power is not allowed at a node without all phases                                 | `LoadGenOnAbsentPhase`             |
| a shunt with status on is not allowed at a node without all phases                                               | `ShuntAtReducedPhaseNode`          |
| the phases of a fault must be present at its node                                                                | `InvalidShortCircuitPhaseOrType`   |
| no symmetric calculation, state estimation or voltage regulator at a node without all phases                     | `ReducedPhaseNotSupported`         |
| the experimental features must be enabled                                                                        | `ExperimentalFeature`              |

A valid result of `validate_input_data` therefore does not guarantee that the calculation succeeds.

## Supported calculations

| Calculation                    | Reduced-phase network                                                                                         |
|--------------------------------|---------------------------------------------------------------------------------------------------------------|
| asymmetric power flow          | supported, with all calculation methods                                                                       |
| asymmetric short circuit       | supported; a short circuit in a network with one- or two-phase components is always calculated asymmetrically |
| symmetric power flow           | not supported (`ReducedPhaseNotSupported`)                                                                    |
| state estimation               | not supported, symmetric or asymmetric (`ReducedPhaseNotSupported`)                                           |

* **Power flow.** Newton-Raphson with its default (linear) initialization is the most robust method.
  Networks with connectors of a near-zero impedance need a different tolerance or a merge; see
  [Networks with near-zero impedance branches](./calculations.md#networks-with-near-zero-impedance-branches).
  The iterative current method handles a phase-to-phase secondary that is not earthed (for example an `ab` → `ab`
  transformer) when a constant-impedance (`const_y`) load is connected to it, because such loads are part of its
  matrix in a network with an absent phase. With only constant-power or constant-current loads there, the
  common-mode voltage of the winding is held only by a very small admittance to earth: the method stays finite but
  does not converge and raises `IterationDiverge` at `max_iterations`. Use Newton-Raphson there.
  Newton-Raphson then gives a common-mode voltage that the loads alone set; with constant-current loads that value may
  be meaningless, because the two loads in series would have to carry different currents.
  The validation warns about such a winding (see [Validation](#validation)).
  The iterative current method always starts from the average-source voltage and ignores `calculation_initialization`.
  The start values of the flat and average-source initializations put leg b of a split-phase node opposite to leg a;
  a transformer that is fed *from* a split-phase node still gets a start that assumes 120° between the phases.
* **Short circuit.** See [Fault](components.md#fault) for the faults that are valid at a node without all phases.
  The default phases of a fault type are not adapted to the node: at a node without them, set `fault_phase`
  explicitly.
  The voltage factor $c$ follows the `u_rated` of the source node, not of the fault location.
  A network whose only special components are delta loads may still use the symmetric short circuit, because loads
  are ignored in short circuits.
* **Not supported:**
  * sources at a node without all phases: a node with a source always has all three phases;
  * voltage regulators at a node without all phases, and on a delta load or generator
    (`UnsupportedLoadGenTypeForVoltageRegulator`).

## C API compatibility

This feature adds attributes to existing components: `phases` to `asym_line`, and `connection` to `sym_load`,
`sym_gen`, `asym_load` and `asym_gen`.
`PGM_create_buffer` returns uninitialised memory, so a C API caller that only sets the attributes it knows leaves the
new attributes as garbage:

* a garbage `connection` of exactly 1 (`delta`) makes an `asym_load` or `asym_gen` a delta load or generator, which
  raises `ExperimentalFeature` in a calculation unless the experimental features are enabled (in which case nothing is
  reported), and makes a `sym_load` or `sym_gen` fail the construction with `InvalidPhaseConnection`;
* a garbage `phases` of 1 to 6 makes the line a one- or two-phase line, which usually fails the construction with
  `InvalidPhaseConnection`, because the entries of the absent phases are provided;
* any other garbage `connection` or `phases` value fails the construction with `MissingCaseForEnumError`; a garbage 0
  is harmless, because it is `wye` for `connection` and `abc` for `phases`.

Call `PGM_buffer_set_nan` on every input buffer before setting its attributes, as the C examples do.
Not provided (the `na` value that `PGM_buffer_set_nan` sets) means `abc` for `phases` and `wye` for `connection`, so
the behaviour of existing models does not change.
See [Backwards compatibility](../advanced_documentation/c-api.md#backwards-compatibility) of the C API.

## Limitations

* The features above are experimental (see [Experimental feature](#experimental-feature)).
* State estimation, symmetric power flow, and sources and voltage regulators at nodes without all phases are not
  supported.
* The automatic tap changer of a `single_phase_transformer` works in asymmetric power flow only, like the transformer
  itself.
* A node has no neutral conductor: a wye load at a split-phase node is connected between its leg and true earth, and
  the two windings of an `a` → `a` transformer cannot share one earth electrode.
* A delta load or generator gives per-phase output only, not per branch; a delta `sym_load` or `sym_gen` is not
  supported.
* The phases of a node are not part of the output.
* The voltage factor $c$ of a short circuit follows the source node, not the fault location.

Besides the small cases for each component, the feature is validated against OpenDSS on two real rural networks with SWER
(SMR8 and KLO14, AusNet Services data from the Team-Nando MV-LV-Networks repository), for power flow and short circuit.
The cases are in `tests/data/power_flow/australian-swer` and `tests/data/short_circuit/australian-swer`, and each has a
README that lists its modelling choices.

The [design of reduced-phase components](../advanced_documentation/multi-phase-and-swer-design.md#later-work) lists the
work planned for later, and also records the decisions behind this feature.
