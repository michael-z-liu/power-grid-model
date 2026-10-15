# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

from copy import copy

import numpy as np
import pytest

from power_grid_model import PowerGridModel
from power_grid_model._core.dataset_definitions import AttributeType as AT, ComponentType as CT, DatasetType
from power_grid_model._core.error_handling import _interpret_error
from power_grid_model._core.power_grid_meta import initialize_array
from power_grid_model.enum import (
    AngleMeasurementType,
    CalculationMethod,
    LoadGenType,
    MeasuredTerminalType,
    PhaseConnection,
    TapChangingStrategy,
    WindingType,
    _ExperimentalFeatures,
)
from power_grid_model.errors import (
    AutomaticTapInputError,
    ConflictID,
    ConflictingAngleMeasurementType,
    ConflictVoltage,
    DanglingConductor,
    ExperimentalFeature,
    IDNotFound,
    IDWrongType,
    InvalidBranch,
    InvalidBranch3,
    InvalidCalculationMethod,
    InvalidMeasuredObject,
    InvalidPhaseConnection,
    InvalidRegulatedObject,
    InvalidTransformerClock,
    IterationDiverge,
    LoadGenOnAbsentPhase,
    MissingCaseForEnumError,
    MissingPhaseForWinding,
    NotObservableError,
    PowerGridBatchError,
    PowerGridError,
    ReducedPhaseNotSupported,
    ShuntAtReducedPhaseNode,
)

from .utils import PowerGridModelWithExt


def test_empty_model():
    model = PowerGridModel.__new__(PowerGridModel)
    with pytest.raises(TypeError):
        model.calculate_power_flow()
    with pytest.raises(TypeError):
        model.copy()
    with pytest.raises(TypeError):
        _n = model.all_component_count
    with pytest.raises(TypeError):
        copy(model)


def test_unknown_component_types():
    model = PowerGridModel(input_data={})
    with pytest.raises(KeyError, match=r"fake_type"):
        model.calculate_power_flow(output_component_types={"fake_type"})


def test_handle_conflict_voltage_error():
    node_input = initialize_array(DatasetType.input, CT.node, 2)
    node_input[AT.id] = [0, 1]
    node_input[AT.u_rated] = [10.0e3, 20.0e3]

    line_input = initialize_array(DatasetType.input, CT.line, 1)
    line_input[AT.id] = [2]
    line_input[AT.from_node] = [0]
    line_input[AT.to_node] = [1]

    with pytest.raises(ConflictVoltage):
        PowerGridModel(input_data={CT.node: node_input, CT.line: line_input})


def test_handle_missing_case_for_enum_error():
    node_input = initialize_array(DatasetType.input, CT.node, 2)
    node_input[AT.id] = [0, 1]
    node_input[AT.u_rated] = [1e4, 4e2]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [2]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [10.0e3]

    transformer_input = initialize_array(DatasetType.input, CT.transformer, 1)
    transformer_input[AT.id] = [3]
    transformer_input[AT.from_node] = [0]
    transformer_input[AT.to_node] = [1]
    transformer_input[AT.from_status] = [1]
    transformer_input[AT.to_status] = [1]
    transformer_input[AT.winding_from] = [2]
    transformer_input[AT.winding_to] = [1]
    transformer_input[AT.clock] = [5]
    transformer_input[AT.tap_side] = [0]

    transformer_tap_regulator_input = initialize_array(DatasetType.input, CT.transformer_tap_regulator, 1)
    transformer_tap_regulator_input[AT.id] = [4]
    transformer_tap_regulator_input[AT.regulated_object] = [3]
    transformer_tap_regulator_input[AT.status] = [1]
    transformer_tap_regulator_input[AT.control_side] = [127]

    with pytest.raises(MissingCaseForEnumError):
        PowerGridModelWithExt(
            input_data={
                CT.node: node_input,
                CT.transformer: transformer_input,
                CT.source: source_input,
                CT.transformer_tap_regulator: transformer_tap_regulator_input,
            }
        )


def test_handle_branch_into_itself_is_valid():
    node_input = initialize_array(DatasetType.input, CT.node, 1)
    node_input[AT.id] = [0]
    node_input[AT.u_rated] = [0.0]

    line_input = initialize_array(DatasetType.input, CT.line, 1)
    line_input[AT.id] = [1]
    line_input[AT.from_node] = [0]
    line_input[AT.to_node] = [0]

    PowerGridModel(input_data={CT.node: node_input, CT.line: line_input})


def test_handle_branch3_into_itself_is_valid():
    node_input = initialize_array(DatasetType.input, CT.node, 1)
    node_input[AT.id] = [0]
    node_input[AT.u_rated] = [0.0]

    three_winding_transformer_input = initialize_array(DatasetType.input, CT.three_winding_transformer, 1)
    three_winding_transformer_input[AT.id] = [1]
    three_winding_transformer_input[AT.node_1] = [0]
    three_winding_transformer_input[AT.node_2] = [0]
    three_winding_transformer_input[AT.node_3] = [0]

    PowerGridModel(
        input_data={
            CT.node: node_input,
            CT.three_winding_transformer: three_winding_transformer_input,
        }
    )


def test_handle_invalid_transformer_clock_error():
    node_input = initialize_array(DatasetType.input, CT.node, 2)
    node_input[AT.id] = [0, 1]
    node_input[AT.u_rated] = [0.0, 0.0]

    transformer_input = initialize_array(DatasetType.input, CT.transformer, 1)
    transformer_input[AT.id] = [2]
    transformer_input[AT.from_node] = [0]
    transformer_input[AT.to_node] = [1]
    transformer_input[AT.clock] = [-1]

    with pytest.raises(InvalidTransformerClock):
        PowerGridModel(input_data={CT.node: node_input, CT.transformer: transformer_input})


@pytest.mark.skip(reason="TODO")
def test_handle_sparse_matrix_error():
    pass


def test_handle_not_observable_error():
    node_input = initialize_array(DatasetType.input, CT.node, 1)
    node_input[AT.id] = [0]
    node_input[AT.u_rated] = [10.0e3]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [1]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]

    sym_load_input = initialize_array(DatasetType.input, CT.sym_load, 1)
    sym_load_input[AT.id] = [2]
    sym_load_input[AT.node] = [0]
    sym_load_input[AT.status] = [1]
    sym_load_input[AT.type] = [LoadGenType.const_power]

    model = PowerGridModel(
        input_data={
            CT.node: node_input,
            CT.source: source_input,
            CT.sym_load: sym_load_input,
        }
    )
    with pytest.raises(NotObservableError):
        model.calculate_state_estimation(calculation_method=CalculationMethod.iterative_linear)


def test_handle_iteration_diverge_error():
    node_input = initialize_array(DatasetType.input, CT.node, 1)
    node_input[AT.id] = [0]
    node_input[AT.u_rated] = [100.0]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [1]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [1.0]
    source_input[AT.sk] = [1000.0]
    source_input[AT.rx_ratio] = [0.0]

    sym_load_input = initialize_array(DatasetType.input, CT.sym_load, 1)
    sym_load_input[AT.id] = [2]
    sym_load_input[AT.node] = [0]
    sym_load_input[AT.status] = [1]
    sym_load_input[AT.type] = [LoadGenType.const_current]
    sym_load_input[AT.p_specified] = [0.0]
    sym_load_input[AT.q_specified] = [500.0]

    model = PowerGridModel(
        input_data={
            CT.node: node_input,
            CT.source: source_input,
            CT.sym_load: sym_load_input,
        }
    )
    with pytest.raises(IterationDiverge):
        model.calculate_power_flow(max_iterations=1, error_tolerance=1.0e-100)


def test_handle_conflict_id_error():
    node_input = initialize_array(DatasetType.input, CT.node, 2)
    node_input[AT.id] = [0, 0]
    node_input[AT.u_rated] = [0.0, 0.0]

    with pytest.raises(ConflictID):
        PowerGridModel(input_data={CT.node: node_input})


def test_handle_id_not_found_error():
    node_input = initialize_array(DatasetType.input, CT.node, 1)
    node_input[AT.id] = [0]
    node_input[AT.u_rated] = [0.0]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [1]
    source_input[AT.node] = [99]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [0.0]

    with pytest.raises(IDNotFound):
        PowerGridModel(input_data={CT.node: node_input, CT.source: source_input})


@pytest.mark.parametrize("sensor_type", [CT.sym_power_sensor, CT.sym_current_sensor])
def test_handle_invalid_measured_object_error(sensor_type):
    node_input = initialize_array(DatasetType.input, CT.node, 2)
    node_input[AT.id] = [0, 1]
    node_input[AT.u_rated] = [0.0, 0.0]

    link_input = initialize_array(DatasetType.input, CT.link, 1)
    link_input[AT.id] = [2]
    link_input[AT.from_node] = [0]
    link_input[AT.to_node] = [1]

    sensor_input = initialize_array(DatasetType.input, sensor_type, 1)
    sensor_input[AT.id] = [3]
    sensor_input[AT.measured_object] = [2]
    sensor_input[AT.measured_terminal_type] = [MeasuredTerminalType.branch_from]

    with pytest.raises(InvalidMeasuredObject):
        PowerGridModel(input_data={CT.node: node_input, CT.link: link_input, sensor_type: sensor_input})


def test_handle_invalid_regulated_object_error():
    node_input = initialize_array(DatasetType.input, CT.node, 2)
    node_input[AT.id] = [0, 1]
    node_input[AT.u_rated] = [1e4, 4e2]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [2]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [10.0e3]

    transformer_tap_regulator_input = initialize_array(DatasetType.input, CT.transformer_tap_regulator, 1)
    transformer_tap_regulator_input[AT.id] = [3]
    transformer_tap_regulator_input[AT.regulated_object] = [2]

    with pytest.raises(InvalidRegulatedObject):
        PowerGridModel(
            input_data={
                CT.node: node_input,
                CT.source: source_input,
                CT.transformer_tap_regulator: transformer_tap_regulator_input,
            }
        )


def test_handle_id_wrong_type_error():
    node_input = initialize_array(DatasetType.input, CT.node, 1)
    node_input[AT.id] = [0]
    node_input[AT.u_rated] = [0.0]

    sym_power_sensor_input = initialize_array(DatasetType.input, CT.sym_power_sensor, 1)
    sym_power_sensor_input[AT.id] = [1]
    sym_power_sensor_input[AT.measured_object] = [0]
    sym_power_sensor_input[AT.measured_terminal_type] = [MeasuredTerminalType.branch_from]

    with pytest.raises(IDWrongType):
        PowerGridModel(input_data={CT.node: node_input, CT.sym_power_sensor: sym_power_sensor_input})


def test_handle_invalid_calculation_method_error():
    node_input = initialize_array(DatasetType.input, CT.node, 1)
    node_input[AT.id] = [0]
    node_input[AT.u_rated] = [10.0e3]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [1]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [10.0e3]

    sym_load_input = initialize_array(DatasetType.input, CT.sym_load, 1)
    sym_load_input[AT.id] = [2]
    sym_load_input[AT.node] = [0]
    sym_load_input[AT.status] = [1]
    sym_load_input[AT.type] = [LoadGenType.const_power]
    sym_load_input[AT.p_specified] = [10.0e3]
    sym_load_input[AT.q_specified] = [1.0e3]

    model = PowerGridModel(
        input_data={
            CT.node: node_input,
            CT.source: source_input,
            CT.sym_load: sym_load_input,
        }
    )
    with pytest.raises(InvalidCalculationMethod):
        model.calculate_power_flow(calculation_method=CalculationMethod.iec60909)


def test_transformer_tap_regulator_control_side_not_closer_to_source():
    node_input = initialize_array(DatasetType.input, CT.node, 2)
    node_input[AT.id] = [0, 1]
    node_input[AT.u_rated] = [150e3, 10e3]

    transformer_input = initialize_array(DatasetType.input, CT.transformer, 1)
    transformer_input[AT.id] = [2]
    transformer_input[AT.from_node] = [0]
    transformer_input[AT.to_node] = [1]
    transformer_input[AT.from_status] = [1]
    transformer_input[AT.to_status] = [1]
    transformer_input[AT.u1] = [1e4]
    transformer_input[AT.u2] = [4e2]
    transformer_input[AT.sn] = [1e5]
    transformer_input[AT.uk] = [0.1]
    transformer_input[AT.pk] = [1e3]
    transformer_input[AT.i0] = [1.0e-6]
    transformer_input[AT.p0] = [0.1]
    transformer_input[AT.winding_from] = [1]
    transformer_input[AT.winding_to] = [1]
    transformer_input[AT.clock] = [0]
    transformer_input[AT.tap_side] = [1]
    transformer_input[AT.tap_pos] = [0]
    transformer_input[AT.tap_nom] = [0]
    transformer_input[AT.tap_min] = [-1]
    transformer_input[AT.tap_max] = [1]
    transformer_input[AT.tap_size] = [100]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [3]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [1.0]

    transformer_tap_regulator_input = initialize_array(DatasetType.input, CT.transformer_tap_regulator, 1)
    transformer_tap_regulator_input[AT.id] = [4]
    transformer_tap_regulator_input[AT.regulated_object] = [2]
    transformer_tap_regulator_input[AT.status] = [1]
    transformer_tap_regulator_input[AT.u_set] = [10]
    transformer_tap_regulator_input[AT.u_band] = [200]
    transformer_tap_regulator_input[AT.control_side] = [0]

    model = PowerGridModel(
        input_data={
            CT.node: node_input,
            CT.transformer: transformer_input,
            CT.source: source_input,
            CT.transformer_tap_regulator: transformer_tap_regulator_input,
        }
    )
    with pytest.raises(AutomaticTapInputError):
        model.calculate_power_flow(tap_changing_strategy=TapChangingStrategy.min_voltage_tap, decode_error=True)


def test_automatic_tap_changing():
    model = PowerGridModel(input_data={})
    model.calculate_power_flow(tap_changing_strategy=TapChangingStrategy.min_voltage_tap)


def test_conflicting_angle_measurement_type() -> None:
    node_input = initialize_array(DatasetType.input, CT.node, 2)
    node_input[AT.id] = [0, 1]
    node_input[AT.u_rated] = [10e3, 10e3]

    line_input = initialize_array(DatasetType.input, CT.line, 1)
    line_input[AT.id] = [2]
    line_input[AT.from_node] = [0]
    line_input[AT.to_node] = [1]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [3]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [1.0]

    sym_voltage_sensor_input = initialize_array(DatasetType.input, CT.sym_voltage_sensor, 1)
    sym_voltage_sensor_input[AT.id] = [4]
    sym_voltage_sensor_input[AT.measured_object] = [0]
    sym_voltage_sensor_input[AT.u_sigma] = [1.0]
    sym_voltage_sensor_input[AT.u_measured] = [1.0]
    sym_voltage_sensor_input[AT.u_angle_measured] = [0.0]

    sym_current_sensor_input = initialize_array(DatasetType.input, CT.sym_current_sensor, 2)
    sym_current_sensor_input[AT.id] = [5, 6]
    sym_current_sensor_input[AT.measured_object] = [2, 2]
    sym_current_sensor_input[AT.measured_terminal_type] = [
        MeasuredTerminalType.branch_from,
        MeasuredTerminalType.branch_from,
    ]
    sym_current_sensor_input[AT.angle_measurement_type] = [
        AngleMeasurementType.local_angle,
        AngleMeasurementType.global_angle,
    ]
    sym_current_sensor_input[AT.i_sigma] = [1.0, 1.0]
    sym_current_sensor_input[AT.i_angle_sigma] = [1.0, 1.0]
    sym_current_sensor_input[AT.i_measured] = [0.0, 0.0]
    sym_current_sensor_input[AT.i_angle_measured] = [0.0, 0.0]

    model = PowerGridModel(
        input_data={
            CT.node: node_input,
            CT.line: line_input,
            CT.source: source_input,
            CT.sym_voltage_sensor: sym_voltage_sensor_input,
            CT.sym_current_sensor: sym_current_sensor_input,
        }
    )

    with pytest.raises(ConflictingAngleMeasurementType):
        model._calculate_state_estimation(decode_error=True)


def test_global_current_measurement_without_voltage_angle() -> None:
    node_input = initialize_array(DatasetType.input, CT.node, 2)
    node_input[AT.id] = [0, 1]
    node_input[AT.u_rated] = [10e3, 10e3]

    line_input = initialize_array(DatasetType.input, CT.line, 1)
    line_input[AT.id] = [2]
    line_input[AT.from_node] = [0]
    line_input[AT.to_node] = [1]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [3]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [1.0]

    sym_voltage_sensor_input = initialize_array(DatasetType.input, CT.sym_voltage_sensor, 1)
    sym_voltage_sensor_input[AT.id] = [4]
    sym_voltage_sensor_input[AT.measured_object] = [0]
    sym_voltage_sensor_input[AT.u_sigma] = [1.0]
    sym_voltage_sensor_input[AT.u_measured] = [0.0]

    sym_current_sensor_input = initialize_array(DatasetType.input, CT.sym_current_sensor, 1)
    sym_current_sensor_input[AT.id] = [5]
    sym_current_sensor_input[AT.measured_object] = [2]
    sym_current_sensor_input[AT.measured_terminal_type] = [
        MeasuredTerminalType.branch_from,
    ]
    sym_current_sensor_input[AT.angle_measurement_type] = [
        AngleMeasurementType.global_angle,
    ]
    sym_current_sensor_input[AT.i_sigma] = [1.0]
    sym_current_sensor_input[AT.i_measured] = [0.0]

    model = PowerGridModel(
        input_data={
            CT.node: node_input,
            CT.line: line_input,
            CT.source: source_input,
            CT.sym_voltage_sensor: sym_voltage_sensor_input,
            CT.sym_current_sensor: sym_current_sensor_input,
        }
    )

    with pytest.raises(NotObservableError):
        model._calculate_state_estimation(decode_error=True)


@pytest.mark.skip(reason="TODO")
def test_handle_power_grid_dataset_error():
    pass


@pytest.mark.skip(reason="TODO")
def test_handle_power_grid_unreachable_error():
    pass


def test_deprecated_invalid_branch_error():
    with pytest.deprecated_call():
        InvalidBranch()


def test_deprecated_invalid_branch3_error():
    with pytest.deprecated_call():
        InvalidBranch3()


@pytest.mark.parametrize(
    ("message", "expected_type"),
    [
        pytest.param(
            "State estimation is not supported in a network with one- or two-phase nodes.\n",
            ReducedPhaseNotSupported,
            id="reduced_phase_state_estimation",
        ),
        pytest.param(
            "Symmetric power flow is not supported in a network with one- or two-phase nodes.\n",
            ReducedPhaseNotSupported,
            id="reduced_phase_symmetric_power_flow",
        ),
        pytest.param(
            "Voltage regulation at a node without all phases is not supported in a network with one- or two-phase "
            "nodes.\n",
            ReducedPhaseNotSupported,
            id="reduced_phase_voltage_regulation",
        ),
        pytest.param(
            "Conductors of branch #2 on phases (bc) are not present at any of its connected nodes (from node phases: "
            "a, to node phases: open). Every conductor must be present at the nodes it connects.\n",
            DanglingConductor,
            id="dangling_conductor_index",
        ),
        pytest.param(
            "Conductors of asym_line with ID 42 on phases (b) are not present at any of its connected nodes (from node "
            "phases: a, to node phases: open). Every conductor must be present at the nodes it connects.\n",
            DanglingConductor,
            id="dangling_conductor_id",
        ),
        pytest.param(
            "Winding of branch #3 has phases (abc) that are not all present at its node (node phases: a). A connected "
            "winding must have all of its phases present.\n",
            MissingPhaseForWinding,
            id="missing_phase_for_winding_index",
        ),
        pytest.param(
            "Winding of branch3 with ID 7 has phases (abc) that are not all present at its node (node phases: a). A "
            "connected winding must have all of its phases present.\n",
            MissingPhaseForWinding,
            id="missing_phase_for_winding_id",
        ),
        pytest.param(
            "Load or generator 12 has a non-zero specified power on phases (bc) that are not present at its node (node "
            "phases: a). Use an asym_load or asym_gen with zero power on the absent phases; a sym_load or sym_gen "
            "puts power on all three phases.\n",
            LoadGenOnAbsentPhase,
            id="load_gen_on_absent_phase",
        ),
        pytest.param(
            "Shunt 500 is connected to a node without all phases. A shunt at a node without all phases is not "
            "supported.\n",
            ShuntAtReducedPhaseNode,
            id="shunt_at_reduced_phase_node",
        ),
    ],
)
def test_interpret_reduced_phase_errors(message: str, expected_type: type[PowerGridError]):
    assert type(_interpret_error(message)) is expected_type


def _reduced_phase_input(line_phases: list[int], line_nodes: list[tuple[int, int]], n_node: int):
    """node 0 is a three-phase source node; each line connects the given nodes with the given phases (a = 1, b = 2)"""
    node_input = initialize_array(DatasetType.input, CT.node, n_node)
    node_input[AT.id] = list(range(n_node))
    node_input[AT.u_rated] = 10e3

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [100]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [1.0]

    line_input = initialize_array(DatasetType.input, CT.asym_line, len(line_phases))
    line_input[AT.id] = [200 + k for k in range(len(line_phases))]
    line_input[AT.from_node] = [nodes[0] for nodes in line_nodes]
    line_input[AT.to_node] = [nodes[1] for nodes in line_nodes]
    line_input[AT.from_status] = 1
    line_input[AT.to_status] = 1
    line_input["phases"] = line_phases
    for k, phases in enumerate(line_phases):
        phase = {PhaseConnection.a: "a", PhaseConnection.b: "b"}[PhaseConnection(phases)]
        line_input[f"r_{phase}{phase}"][k] = 0.5
        line_input[f"x_{phase}{phase}"][k] = 0.4
        line_input[f"c_{phase}{phase}"][k] = 1e-8

    load_input = initialize_array(DatasetType.input, CT.asym_load, 1)
    load_input[AT.id] = [300]
    load_input[AT.node] = [1]
    load_input[AT.status] = [1]
    load_input[AT.type] = [LoadGenType.const_power]
    load_input[AT.p_specified] = [[1e3, 0.0, 0.0]]
    load_input[AT.q_specified] = [[0.0, 0.0, 0.0]]

    return {CT.node: node_input, CT.source: source_input, CT.asym_line: line_input, CT.asym_load: load_input}


def test_reduced_phase_calculations():
    model = PowerGridModelWithExt(_reduced_phase_input([PhaseConnection.a], [(0, 1)], 2))
    enabled = _ExperimentalFeatures.enabled

    result = model.calculate_power_flow_with_ext(symmetric=False, experimental_features=enabled)
    assert result[CT.node][AT.u][1][0] > 0.0
    assert result[CT.node][AT.u][1][1] == 0.0
    assert result[CT.node][AT.u][1][2] == 0.0

    with pytest.raises(ExperimentalFeature):
        model.calculate_power_flow(symmetric=False)
    with pytest.raises(ReducedPhaseNotSupported, match="Symmetric power flow"):
        model.calculate_power_flow_with_ext(symmetric=True, experimental_features=enabled)
    for symmetric in (True, False):
        with pytest.raises(ReducedPhaseNotSupported, match="State estimation"):
            model.calculate_state_estimation_with_ext(symmetric=symmetric, experimental_features=enabled)

    # a short circuit without faults runs through the asymmetric calculation
    sc_result = model.calculate_short_circuit_with_ext(experimental_features=enabled)
    assert sc_result[CT.node][AT.u][1][0] > 0.0
    assert sc_result[CT.node][AT.u][1][1] == 0.0


def test_reduced_phase_dangling_conductor():
    # a b line off an a node
    model = PowerGridModelWithExt(_reduced_phase_input([PhaseConnection.a, PhaseConnection.b], [(0, 1), (1, 2)], 3))
    with pytest.raises(DanglingConductor, match="branch with ID 201"):
        model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)


def test_reduced_phase_missing_phase_for_winding():
    input_data = _reduced_phase_input([PhaseConnection.a], [(0, 1)], 3)
    input_data[CT.node][AT.u_rated][2] = 400.0

    transformer_input = initialize_array(DatasetType.input, CT.transformer, 1)
    transformer_input[AT.id] = [400]
    transformer_input[AT.from_node] = [1]
    transformer_input[AT.to_node] = [2]
    transformer_input[AT.from_status] = [1]
    transformer_input[AT.to_status] = [1]
    transformer_input[AT.u1] = [10e3]
    transformer_input[AT.u2] = [400.0]
    transformer_input[AT.sn] = [1e5]
    transformer_input[AT.uk] = [0.05]
    transformer_input[AT.pk] = [1e3]
    transformer_input[AT.i0] = [0.01]
    transformer_input[AT.p0] = [100.0]
    transformer_input[AT.winding_from] = [WindingType.wye_n]
    transformer_input[AT.winding_to] = [WindingType.wye_n]
    transformer_input[AT.clock] = [0]
    transformer_input[AT.tap_side] = [0]
    transformer_input[AT.tap_pos] = [0]
    transformer_input[AT.tap_min] = [0]
    transformer_input[AT.tap_max] = [0]
    transformer_input[AT.tap_size] = [100.0]
    input_data[CT.transformer] = transformer_input

    model = PowerGridModelWithExt(input_data)
    with pytest.raises(MissingPhaseForWinding, match="branch with ID 400"):
        model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)


def test_reduced_phase_missing_phase_for_three_winding_transformer():
    input_data = _reduced_phase_input([PhaseConnection.a], [(0, 1)], 4)
    input_data[CT.node][AT.u_rated][2:] = 400.0

    transformer_input = initialize_array(DatasetType.input, CT.three_winding_transformer, 1)
    transformer_input[AT.id] = [400]
    transformer_input[AT.node_1] = [1]
    transformer_input[AT.node_2] = [2]
    transformer_input[AT.node_3] = [3]
    transformer_input[AT.status_1] = [1]
    transformer_input[AT.status_2] = [1]
    transformer_input[AT.status_3] = [1]
    transformer_input[AT.u1] = [10e3]
    transformer_input[AT.u2] = [400.0]
    transformer_input[AT.u3] = [400.0]
    transformer_input[AT.sn_1] = [1e5]
    transformer_input[AT.sn_2] = [1e5]
    transformer_input[AT.sn_3] = [1e5]
    transformer_input[AT.uk_12] = [0.05]
    transformer_input[AT.uk_13] = [0.05]
    transformer_input[AT.uk_23] = [0.05]
    transformer_input[AT.pk_12] = [1e3]
    transformer_input[AT.pk_13] = [1e3]
    transformer_input[AT.pk_23] = [1e3]
    transformer_input[AT.i0] = [0.01]
    transformer_input[AT.p0] = [100.0]
    transformer_input[AT.winding_1] = [WindingType.wye_n]
    transformer_input[AT.winding_2] = [WindingType.wye_n]
    transformer_input[AT.winding_3] = [WindingType.wye_n]
    transformer_input[AT.clock_12] = [0]
    transformer_input[AT.clock_13] = [0]
    transformer_input[AT.tap_side] = [0]
    transformer_input[AT.tap_pos] = [0]
    transformer_input[AT.tap_min] = [0]
    transformer_input[AT.tap_max] = [0]
    transformer_input[AT.tap_size] = [100.0]
    input_data[CT.three_winding_transformer] = transformer_input

    model = PowerGridModelWithExt(input_data)
    with pytest.raises(MissingPhaseForWinding, match="branch3 with ID 400"):
        model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)


def test_reduced_phase_invalid_phase_connection():
    input_data = _reduced_phase_input([PhaseConnection.a], [(0, 1)], 2)
    input_data[CT.asym_line]["r_bb"] = [0.5]
    with pytest.raises(InvalidPhaseConnection, match="r_bb"):
        PowerGridModel(input_data)

    input_data = _reduced_phase_input([PhaseConnection.a], [(0, 1)], 2)
    input_data[CT.asym_line]["phases"] = [9]
    with pytest.raises(MissingCaseForEnumError):
        PowerGridModel(input_data)


def _single_phase_node_input(load_type: CT, p: float | list[float], q: float | list[float], status: int = 1):
    """node 0 (abc, source) --a-- node 1 (a) with one load of the given type at node 1"""
    input_data = _reduced_phase_input([PhaseConnection.a], [(0, 1)], 2)
    load_input = initialize_array(DatasetType.input, load_type, 1)
    load_input[AT.id] = [300]
    load_input[AT.node] = [1]
    load_input[AT.status] = [status]
    load_input[AT.type] = [LoadGenType.const_power]
    load_input[AT.p_specified] = [p]
    load_input[AT.q_specified] = [q]
    del input_data[CT.asym_load]
    input_data[load_type] = load_input
    return input_data


@pytest.mark.parametrize(
    ("load_type", "p", "q"),
    [
        pytest.param(CT.sym_load, 0.0, 0.0, id="sym_load_zero"),
        pytest.param(CT.asym_load, [1e3, 0.0, 0.0], [1e2, 0.0, 0.0], id="asym_load_zero_absent"),
        pytest.param(CT.asym_load, [1e3, np.nan, np.nan], [1e2, np.nan, np.nan], id="asym_load_nan_absent"),
        pytest.param(CT.asym_gen, [1e3, 0.0, 0.0], [0.0, 0.0, 0.0], id="asym_gen_zero_absent"),
    ],
)
def test_reduced_phase_load_power_on_present_phases(load_type, p, q):
    model = PowerGridModelWithExt(_single_phase_node_input(load_type, p, q))
    model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)


@pytest.mark.parametrize(
    ("load_type", "p", "q", "status"),
    [
        pytest.param(CT.sym_load, 1e3, 0.0, 1, id="sym_load"),
        pytest.param(CT.sym_gen, 0.0, 1e2, 1, id="sym_gen_q"),
        pytest.param(CT.asym_load, [1e3, 1e2, 0.0], [0.0, 0.0, 0.0], 1, id="asym_load_p_on_b"),
        pytest.param(CT.asym_load, [1e3, 0.0, 0.0], [0.0, 0.0, 1e2], 1, id="asym_load_q_on_c"),
    ],
)
def test_reduced_phase_load_power_on_absent_phases(load_type, p, q, status):
    model = PowerGridModelWithExt(_single_phase_node_input(load_type, p, q, status))
    with pytest.raises(LoadGenOnAbsentPhase, match="Load or generator 300 "):
        model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)


def test_reduced_phase_load_power_on_absent_phases_status_off():
    model = PowerGridModelWithExt(_single_phase_node_input(CT.sym_load, 1e3, 1e2, status=0))
    model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)


def test_reduced_phase_load_power_on_absent_phases_batch():
    model = PowerGridModelWithExt(_single_phase_node_input(CT.sym_load, 0.0, 0.0))

    update = initialize_array(DatasetType.update, CT.sym_load, (2, 1))
    update[AT.id] = 300
    update[AT.p_specified] = [[1e3], [0.0]]
    result = model.calculate_power_flow_with_ext(
        symmetric=False,
        update_data={CT.sym_load: update},
        continue_on_batch_error=True,
        experimental_features=_ExperimentalFeatures.enabled,
    )
    assert result[CT.node][AT.u][1, 1, 0] > 0.0

    with pytest.raises(PowerGridBatchError) as error:
        model.calculate_power_flow_with_ext(
            symmetric=False, update_data={CT.sym_load: update}, experimental_features=_ExperimentalFeatures.enabled
        )
    assert error.value.failed_scenarios.tolist() == [0]
    assert error.value.errors[0].__class__ is LoadGenOnAbsentPhase


def test_reduced_phase_load_power_on_absent_phases_fixed_in_every_scenario():
    # the cache run with the base input fails, but every scenario provides a valid value
    model = PowerGridModelWithExt(_single_phase_node_input(CT.sym_load, 1e3, 0.0))

    update = initialize_array(DatasetType.update, CT.sym_load, (2, 1))
    update[AT.id] = 300
    update[AT.p_specified] = 0.0
    update[AT.status] = [[1], [0]]
    result = model.calculate_power_flow_with_ext(
        symmetric=False, update_data={CT.sym_load: update}, experimental_features=_ExperimentalFeatures.enabled
    )
    assert result[CT.node][AT.u].shape == (2, 2, 3)


@pytest.mark.parametrize("calculation", ["power_flow", "short_circuit"])
def test_reduced_phase_shunt(calculation):
    input_data = _reduced_phase_input([PhaseConnection.a], [(0, 1)], 2)
    shunt_input = initialize_array(DatasetType.input, CT.shunt, 1)
    shunt_input[AT.id] = [500]
    shunt_input[AT.node] = [1]
    shunt_input[AT.status] = [1]
    shunt_input[AT.g1] = [0.0]
    shunt_input[AT.b1] = [1e-4]
    shunt_input[AT.g0] = [0.0]
    shunt_input[AT.b0] = [1e-4]
    input_data[CT.shunt] = shunt_input

    def calculate(model):
        if calculation == "power_flow":
            return model.calculate_power_flow_with_ext(
                symmetric=False, experimental_features=_ExperimentalFeatures.enabled
            )
        return model.calculate_short_circuit_with_ext(experimental_features=_ExperimentalFeatures.enabled)

    with pytest.raises(ShuntAtReducedPhaseNode, match="Shunt 500 is connected to a node without all phases"):
        calculate(PowerGridModelWithExt(input_data))

    input_data[CT.shunt][AT.status] = [0]
    calculate(PowerGridModelWithExt(input_data))


def _single_phase_node_shunt_input(status: int = 1):
    """node 0 (abc, source) --a-- node 1 (a) with a shunt at node 1"""
    input_data = _reduced_phase_input([PhaseConnection.a], [(0, 1)], 2)
    shunt_input = initialize_array(DatasetType.input, CT.shunt, 1)
    shunt_input[AT.id] = [500]
    shunt_input[AT.node] = [1]
    shunt_input[AT.status] = [status]
    shunt_input[AT.g1] = [0.0]
    shunt_input[AT.b1] = [1e-4]
    shunt_input[AT.g0] = [0.0]
    shunt_input[AT.b0] = [1e-4]
    input_data[CT.shunt] = shunt_input
    return input_data


def test_reduced_phase_shunt_switched_off_in_every_scenario():
    # the cache run with the base input sees an energised shunt, but every scenario switches it off
    model = PowerGridModelWithExt(_single_phase_node_shunt_input())

    update = initialize_array(DatasetType.update, CT.shunt, (2, 1))
    update[AT.id] = 500
    update[AT.status] = 0
    result = model.calculate_power_flow_with_ext(
        symmetric=False, update_data={CT.shunt: update}, experimental_features=_ExperimentalFeatures.enabled
    )
    assert result[CT.node][AT.u].shape == (2, 2, 3)

    # a single, non-batch calculation on the same model still raises
    with pytest.raises(ShuntAtReducedPhaseNode, match="Shunt 500 is connected to a node without all phases"):
        model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)


def test_reduced_phase_shunt_left_on_in_one_scenario():
    model = PowerGridModelWithExt(_single_phase_node_shunt_input())

    update = initialize_array(DatasetType.update, CT.shunt, (2, 1))
    update[AT.id] = 500
    update[AT.status] = [[1], [0]]
    result = model.calculate_power_flow_with_ext(
        symmetric=False,
        update_data={CT.shunt: update},
        continue_on_batch_error=True,
        experimental_features=_ExperimentalFeatures.enabled,
    )
    assert result[CT.node][AT.u][1, 1, 0] > 0.0
    assert model.batch_error is not None
    assert model.batch_error.failed_scenarios.tolist() == [0]
    assert model.batch_error.errors[0].__class__ is ShuntAtReducedPhaseNode

    with pytest.raises(PowerGridBatchError) as error:
        model.calculate_power_flow_with_ext(
            symmetric=False, update_data={CT.shunt: update}, experimental_features=_ExperimentalFeatures.enabled
        )
    assert error.value.failed_scenarios.tolist() == [0]
    assert error.value.errors[0].__class__ is ShuntAtReducedPhaseNode
