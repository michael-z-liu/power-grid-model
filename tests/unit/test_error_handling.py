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
    BranchSide,
    CalculationInitialization,
    CalculationMethod,
    FaultPhase,
    FaultType,
    LoadGenConnection,
    LoadGenType,
    MeasuredTerminalType,
    PhaseConnection,
    TapChangingStrategy,
    WindingConnection1P,
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
    InvalidShortCircuitPhaseOrType,
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
    UnsupportedLoadGenTypeForVoltageRegulator,
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


def _regulated_transformer_behind_branch(branch_type: CT):
    """source -> branch_type -> regulated 10/0.4 kV transformer (u_set 420 V, u_band 8 V) -> 28 kW load"""
    node_input = initialize_array(DatasetType.input, CT.node, 3)
    node_input[AT.id] = [1, 2, 3]
    node_input[AT.u_rated] = [10e3, 10e3, 400.0]

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [4]
    source_input[AT.node] = [1]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [1.0]

    branch_input = initialize_array(DatasetType.input, branch_type, 1)
    branch_input[AT.id] = [5]
    branch_input[AT.from_node] = [1]
    branch_input[AT.to_node] = [2]
    branch_input[AT.from_status] = [1]
    branch_input[AT.to_status] = [1]
    if branch_type == CT.line:
        branch_input[AT.r1] = branch_input[AT.r0] = [0.1]
        branch_input[AT.x1] = branch_input[AT.x0] = [0.1]
        branch_input[AT.c1] = branch_input[AT.c0] = [0.0]
        branch_input[AT.tan1] = branch_input[AT.tan0] = [0.0]
    elif branch_type == CT.asym_line:
        for attribute in (AT.r_aa, AT.r_bb, AT.r_cc, AT.x_aa, AT.x_bb, AT.x_cc):
            branch_input[attribute] = [0.1]
        for attribute in (AT.r_ba, AT.r_ca, AT.r_cb, AT.x_ba, AT.x_ca, AT.x_cb, AT.c0, AT.c1):
            branch_input[attribute] = [0.0]
    else:
        branch_input[AT.r1] = [0.1]
        branch_input[AT.x1] = [0.1]
        branch_input[AT.g1] = [0.0]
        branch_input[AT.b1] = [0.0]

    transformer_input = initialize_array(DatasetType.input, CT.transformer, 1)
    transformer_input[AT.id] = [6]
    transformer_input[AT.from_node] = [2]
    transformer_input[AT.to_node] = [3]
    transformer_input[AT.from_status] = [1]
    transformer_input[AT.to_status] = [1]
    transformer_input[AT.u1] = [10e3]
    transformer_input[AT.u2] = [400.0]
    transformer_input[AT.sn] = [100e3]
    transformer_input[AT.uk] = [0.1]
    transformer_input[AT.pk] = [1e3]
    transformer_input[AT.i0] = [1.0e-6]
    transformer_input[AT.p0] = [0.1]
    transformer_input[AT.winding_from] = [WindingType.wye_n]
    transformer_input[AT.winding_to] = [WindingType.wye_n]
    transformer_input[AT.clock] = [0]
    transformer_input[AT.tap_side] = [BranchSide.from_side]
    transformer_input[AT.tap_pos] = [0]
    transformer_input[AT.tap_nom] = [0]
    transformer_input[AT.tap_min] = [-10]
    transformer_input[AT.tap_max] = [10]
    transformer_input[AT.tap_size] = [100.0]

    load_input = initialize_array(DatasetType.input, CT.sym_load, 1)
    load_input[AT.id] = [7]
    load_input[AT.node] = [3]
    load_input[AT.status] = [1]
    load_input[AT.type] = [LoadGenType.const_power]
    load_input[AT.p_specified] = [28e3]
    load_input[AT.q_specified] = [5.6e3]

    regulator_input = initialize_array(DatasetType.input, CT.transformer_tap_regulator, 1)
    regulator_input[AT.id] = [8]
    regulator_input[AT.regulated_object] = [6]
    regulator_input[AT.status] = [1]
    regulator_input[AT.control_side] = [BranchSide.to_side]
    regulator_input[AT.u_set] = [420.0]
    regulator_input[AT.u_band] = [8.0]
    regulator_input[AT.line_drop_compensation_r] = [0.0]
    regulator_input[AT.line_drop_compensation_x] = [0.0]

    return {
        CT.node: node_input,
        CT.source: source_input,
        branch_type: branch_input,
        CT.transformer: transformer_input,
        CT.sym_load: load_input,
        CT.transformer_tap_regulator: regulator_input,
    }


@pytest.mark.parametrize(
    ("branch_type", "symmetric"),
    [
        (CT.line, True),
        (CT.line, False),
        (CT.asym_line, True),
        (CT.asym_line, False),
        (CT.generic_branch, True),
    ],
)
def test_automatic_tap_changing_behind_any_branch(branch_type: CT, symmetric: bool):
    # A regulated transformer reached from the source only through an asym_line or a generic_branch used to be left
    # out of the transformer ranking and silently not regulated (tap_pos na, u 396.5 V outside the band).
    model = PowerGridModel(_regulated_transformer_behind_branch(branch_type))
    result = model.calculate_power_flow(symmetric=symmetric, tap_changing_strategy=TapChangingStrategy.any_valid_tap)

    expected_tap_pos = -5
    u_set = 420.0
    u_band = 8.0
    assert result[CT.transformer_tap_regulator][AT.tap_pos][0] == expected_tap_pos
    u_low_voltage_node = result[CT.node][AT.u][2]
    u_controlled = u_low_voltage_node if symmetric else np.mean(u_low_voltage_node) * np.sqrt(3)
    assert u_controlled == pytest.approx(417.7, abs=0.05)
    assert u_set - u_band / 2 <= u_controlled <= u_set + u_band / 2


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
            "Load or generator 12 has a non-zero specified power on phases (ca) that are not present at its node (node "
            "phases: ab). A delta (phase-to-phase) entry needs both of its phases: ab needs a and b, bc needs b and c, "
            "ca needs c and a.\n",
            LoadGenOnAbsentPhase,
            id="load_gen_on_absent_phase_delta",
        ),
        pytest.param(
            "Symmetric power flow is not supported in a network with delta-connected loads or generators.\n",
            ReducedPhaseNotSupported,
            id="delta_load_gen_symmetric_power_flow",
        ),
        pytest.param(
            "State estimation is not supported in a network with delta-connected loads or generators.\n",
            ReducedPhaseNotSupported,
            id="delta_load_gen_state_estimation",
        ),
        pytest.param(
            "Automatic tap changer has invalid configuration. The single_phase_transformer 152 regulated by "
            "transformer_tap_regulator 402 needs both tap_min and tap_max.",
            AutomaticTapInputError,
            id="single_phase_transformer_regulator_tap_limits",
        ),
        pytest.param(
            "Shunt 500 is connected to a node without all phases. A shunt at a node without all phases is not "
            "supported.\n",
            ShuntAtReducedPhaseNode,
            id="shunt_at_reduced_phase_node",
        ),
        pytest.param(
            "The short circuit phases (b) of fault 601 are not all present at its node 4 (node phases: a); this short "
            "circuit type cannot be applied there\n",
            InvalidShortCircuitPhaseOrType,
            id="short_circuit_phases_absent_at_fault_node",
        ),
        pytest.param(
            "The short circuit phases (a) of fault 601 are not all present at its node 4 (node phases: b); this short "
            "circuit type cannot be applied there; the phases are the default of the fault type, so set fault_phase "
            "explicitly\n",
            InvalidShortCircuitPhaseOrType,
            id="short_circuit_default_phases_absent_at_fault_node",
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


def _missing_phase_for_winding_input():
    """a three-phase transformer 400 at the a node 1 behind the a line 200"""
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
    return input_data


def test_reduced_phase_missing_phase_for_winding():
    model = PowerGridModelWithExt(_missing_phase_for_winding_input())
    with pytest.raises(MissingPhaseForWinding, match="branch with ID 400"):
        model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)


@pytest.mark.parametrize(
    ("input_data", "branch", "branch_id", "error_type"),
    [
        pytest.param(
            _reduced_phase_input([PhaseConnection.a, PhaseConnection.b], [(0, 1), (1, 2)], 3),
            CT.asym_line,
            201,
            DanglingConductor,
            id="dangling_conductor",
        ),
        pytest.param(_missing_phase_for_winding_input(), CT.transformer, 400, MissingPhaseForWinding, id="winding"),
    ],
)
def test_reduced_phase_topology_error_fixed_in_the_update_data(input_data, branch, branch_id, error_type):
    # the cache run with the base input sees the invalid branch; a scenario that switches it off is valid (as for a
    # shunt, D4) and a scenario that leaves it on fails on its own, with the branch ID
    model = PowerGridModelWithExt(input_data)
    enabled = _ExperimentalFeatures.enabled

    update = initialize_array(DatasetType.update, branch, (2, 1))
    update[AT.id] = branch_id
    update[AT.from_status] = 0
    update[AT.to_status] = 0
    result = model.calculate_power_flow_with_ext(
        symmetric=False, update_data={branch: update}, experimental_features=enabled
    )
    assert (result[CT.node][AT.u][:, 1, 0] > 0.0).all()
    assert (result[CT.node][AT.u][:, 1, 1:] == 0.0).all()

    update[AT.from_status] = [[0], [1]]
    update[AT.to_status] = [[0], [1]]
    with pytest.raises(PowerGridBatchError) as error:
        model.calculate_power_flow_with_ext(
            symmetric=False, update_data={branch: update}, experimental_features=enabled
        )
    assert error.value.failed_scenarios.tolist() == [1]
    assert error.value.errors[0].__class__ is error_type
    assert f"branch with ID {branch_id}" in error.value.error_messages[0]

    # a single, non-batch calculation on the same model still raises
    with pytest.raises(error_type, match=f"branch with ID {branch_id}"):
        model.calculate_power_flow_with_ext(symmetric=False, experimental_features=enabled)


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


def _single_phase_transformer_input(behind_line: bool = False):
    """node 0 (abc, 22 kV, source) -- single_phase_transformer 400 ab -> a (22 kV / 12.7 kV, earth electrode 5 ohm)
    -- last node (SWER, u_rated 22 kV) with a load on phase a; behind_line: the from node is the a node 1 behind
    asym_line 200"""
    from_node = 1 if behind_line else 0
    n_node = 3 if behind_line else 2
    input_data = _reduced_phase_input([PhaseConnection.a], [(0, 1)], n_node)
    input_data[CT.node][AT.u_rated] = 22e3
    if not behind_line:
        del input_data[CT.asym_line]
    input_data[CT.asym_load][AT.node] = [n_node - 1]
    input_data[CT.asym_load][AT.p_specified] = [[50e3, 0.0, 0.0]]

    transformer_input = initialize_array(DatasetType.input, CT.single_phase_transformer, 1)
    transformer_input[AT.id] = [400]
    transformer_input[AT.from_node] = [from_node]
    transformer_input[AT.to_node] = [n_node - 1]
    transformer_input[AT.from_status] = [1]
    transformer_input[AT.to_status] = [1]
    transformer_input[AT.u1] = [22e3]
    transformer_input[AT.u2] = [12.7e3]
    transformer_input[AT.sn] = [200e3]
    transformer_input[AT.uk] = [0.04]
    transformer_input[AT.pk] = [2e3]
    transformer_input[AT.i0] = [0.0]
    transformer_input[AT.p0] = [0.0]
    transformer_input[AT.from_connection] = [WindingConnection1P.ab]
    transformer_input[AT.to_connection] = [WindingConnection1P.a]
    transformer_input[AT.r_grounding_to] = [5.0]
    input_data[CT.single_phase_transformer] = transformer_input
    return input_data


def test_winding_connection_1p_exported():
    assert [(member.name, member.value) for member in WindingConnection1P] == [
        ("a", 1),
        ("b", 2),
        ("c", 3),
        ("ab", 4),
        ("ac", 5),
        ("bc", 6),
        ("ab_center_tapped", 7),
    ]


def test_single_phase_transformer_calculations():
    model = PowerGridModelWithExt(_single_phase_transformer_input())
    enabled = _ExperimentalFeatures.enabled

    result = model.calculate_power_flow_with_ext(symmetric=False, experimental_features=enabled)
    u = result[CT.node][AT.u][1]
    assert u[0] > 0.0
    assert u[1] == 0.0
    assert u[2] == 0.0
    # the SWER voltage follows U_a - U_b, which leads U_a by 30 degrees
    assert result[CT.node][AT.u_angle][1][0] == pytest.approx(np.pi / 6.0, abs=0.05)
    p_to = result[CT.single_phase_transformer][AT.p_to][0]
    assert p_to[0] < 0.0
    assert p_to[1] == 0.0
    assert p_to[2] == 0.0
    assert result[CT.single_phase_transformer][AT.loading][0] > 0.0

    with pytest.raises(ExperimentalFeature):
        model.calculate_power_flow(symmetric=False)
    with pytest.raises(ReducedPhaseNotSupported, match="Symmetric power flow"):
        model.calculate_power_flow_with_ext(symmetric=True, experimental_features=enabled)

    fault_input = initialize_array(DatasetType.input, CT.fault, 1)
    fault_input[AT.id] = [600]
    fault_input[AT.status] = [1]
    fault_input[AT.fault_type] = [FaultType.single_phase_to_ground]
    fault_input[AT.fault_phase] = [FaultPhase.a]
    fault_input[AT.fault_object] = [1]
    input_data = _single_phase_transformer_input()
    input_data[CT.fault] = fault_input
    sc_result = PowerGridModelWithExt(input_data).calculate_short_circuit_with_ext(experimental_features=enabled)
    assert sc_result[CT.fault][AT.i_f][0][0] > 0.0
    assert sc_result[CT.single_phase_transformer][AT.i_to][0][0] > 0.0


def test_single_phase_transformer_missing_phase_for_winding():
    # an ab primary at the a node behind asym_line 200
    model = PowerGridModelWithExt(_single_phase_transformer_input(behind_line=True))
    with pytest.raises(MissingPhaseForWinding, match="branch with ID 400"):
        model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)


def test_single_phase_transformer_invalid_input():
    input_data = _single_phase_transformer_input()
    input_data[CT.single_phase_transformer][AT.r_grounding_from] = [5.0]
    with pytest.raises(InvalidPhaseConnection, match="r_grounding_from"):
        PowerGridModel(input_data)

    input_data = _single_phase_transformer_input()
    input_data[CT.single_phase_transformer][AT.from_connection] = [0]
    with pytest.raises(MissingCaseForEnumError):
        PowerGridModel(input_data)

    input_data = _single_phase_transformer_input()
    input_data[CT.single_phase_transformer][AT.clock] = [3]
    with pytest.raises(InvalidTransformerClock):
        PowerGridModel(input_data)


def test_single_phase_transformer_tap_batch():
    # the tap limits are not provided, so the updated tap positions are not limited
    input_data = _single_phase_transformer_input()
    input_data[CT.single_phase_transformer][AT.tap_size] = [220.0]
    model = PowerGridModelWithExt(input_data)

    update = initialize_array(DatasetType.update, CT.single_phase_transformer, (3, 1))
    update[AT.id] = 400
    update[AT.tap_pos] = [[-2], [0], [2]]
    result = model.calculate_power_flow_with_ext(
        symmetric=False,
        update_data={CT.single_phase_transformer: update},
        experimental_features=_ExperimentalFeatures.enabled,
    )
    u = result[CT.node][AT.u][:, 1, 0]
    # a higher tap position on the from side raises u1 and lowers the secondary voltage in proportion
    assert u[0] / u[1] == pytest.approx(22e3 / 21.56e3, rel=1e-3)
    assert u[1] / u[2] == pytest.approx(22.44e3 / 22e3, rel=1e-3)


def _regulated_single_phase_transformer_input(tap_max_provided: bool = True, regulator_status: int = 1):
    """_single_phase_transformer_input with from-side taps -4..4 of 275 V and transformer_tap_regulator 500 controlling
    the a winding at the to side (u_set 12975 V, u_band 80 V)"""
    input_data = _single_phase_transformer_input()
    transformer_input = input_data[CT.single_phase_transformer]
    transformer_input[AT.tap_side] = [BranchSide.from_side]
    transformer_input[AT.tap_pos] = [0]
    transformer_input[AT.tap_nom] = [0]
    transformer_input[AT.tap_min] = [-4]
    if tap_max_provided:
        transformer_input[AT.tap_max] = [4]
    transformer_input[AT.tap_size] = [275.0]

    regulator_input = initialize_array(DatasetType.input, CT.transformer_tap_regulator, 1)
    regulator_input[AT.id] = [500]
    regulator_input[AT.regulated_object] = [400]
    regulator_input[AT.status] = [regulator_status]
    regulator_input[AT.control_side] = [BranchSide.to_side]
    regulator_input[AT.u_set] = [12975.0]
    regulator_input[AT.u_band] = [80.0]
    input_data[CT.transformer_tap_regulator] = regulator_input
    return input_data


@pytest.mark.parametrize(
    "strategy",
    [
        TapChangingStrategy.any_valid_tap,
        TapChangingStrategy.fast_any_tap,
        TapChangingStrategy.min_voltage_tap,
        TapChangingStrategy.max_voltage_tap,
    ],
)
def test_single_phase_transformer_regulated(strategy: TapChangingStrategy):
    # fixed-tap scan of |U_a| at the SWER node: 12809.2 V at tap -1, 12974.7 V at tap -2, 13144.6 V at tap -3
    model = PowerGridModelWithExt(_regulated_single_phase_transformer_input())
    result = model.calculate_power_flow_with_ext(
        symmetric=False, tap_changing_strategy=strategy, experimental_features=_ExperimentalFeatures.enabled
    )
    expected_tap_pos = -2
    u_set = 12975.0
    u_band = 80.0
    assert result[CT.transformer_tap_regulator][AT.tap_pos][0] == expected_tap_pos
    u_a = result[CT.node][AT.u][1][0]
    assert u_a == pytest.approx(12974.7, abs=0.05)
    assert u_set - u_band / 2 <= u_a <= u_set + u_band / 2

    with pytest.raises(ExperimentalFeature):
        model.calculate_power_flow(symmetric=False, tap_changing_strategy=strategy)
    with pytest.raises(ReducedPhaseNotSupported, match="Symmetric power flow"):
        model.calculate_power_flow_with_ext(
            symmetric=True, tap_changing_strategy=strategy, experimental_features=_ExperimentalFeatures.enabled
        )


@pytest.mark.parametrize("regulator_status", [1, 0])
def test_single_phase_transformer_regulator_needs_tap_limits(regulator_status: int):
    input_data = _regulated_single_phase_transformer_input(tap_max_provided=False, regulator_status=regulator_status)
    with pytest.raises(
        AutomaticTapInputError,
        match=r"The single_phase_transformer 400 regulated by transformer_tap_regulator 500 needs both tap_min and "
        r"tap_max\.",
    ):
        PowerGridModel(input_data)


def test_single_phase_transformer_status_batch():
    # a batch that toggles the status of each side gives the same results as separate calculations
    statuses = [(1, 1), (0, 1), (1, 0), (0, 0), (1, 1)]
    enabled = _ExperimentalFeatures.enabled
    input_data = _single_phase_transformer_input()
    # with a magnetising admittance the connected winding of a half-open transformer sees the open-end shunt
    input_data[CT.single_phase_transformer][AT.i0] = [0.01]
    input_data[CT.single_phase_transformer][AT.p0] = [200.0]

    update = initialize_array(DatasetType.update, CT.single_phase_transformer, (len(statuses), 1))
    update[AT.id] = 400
    update[AT.from_status] = [[from_status] for from_status, _ in statuses]
    update[AT.to_status] = [[to_status] for _, to_status in statuses]
    batch_result = PowerGridModelWithExt(input_data).calculate_power_flow_with_ext(
        symmetric=False,
        update_data={CT.single_phase_transformer: update},
        experimental_features=enabled,
    )

    for scenario, (from_status, to_status) in enumerate(statuses):
        single_input = {component: data.copy() for component, data in input_data.items()}
        single_input[CT.single_phase_transformer][AT.from_status] = [from_status]
        single_input[CT.single_phase_transformer][AT.to_status] = [to_status]
        single_result = PowerGridModelWithExt(single_input).calculate_power_flow_with_ext(
            symmetric=False, experimental_features=enabled
        )
        for component, attributes in (
            (CT.node, (AT.u, AT.u_angle, AT.p, AT.q)),
            (CT.single_phase_transformer, (AT.p_from, AT.q_from, AT.i_from, AT.p_to, AT.q_to, AT.i_to, AT.loading)),
            (CT.asym_load, (AT.p, AT.q, AT.i)),
        ):
            for attribute in attributes:
                np.testing.assert_allclose(
                    batch_result[component][attribute][scenario],
                    single_result[component][attribute],
                    rtol=1e-12,
                    atol=1e-9,
                    err_msg=f"scenario {scenario} {component} {attribute}",
                )
    # an open side carries no current and the energized side of a half-open transformer only the open-end shunt
    i_from = batch_result[CT.single_phase_transformer][AT.i_from][:, 0, 0]
    i_to = batch_result[CT.single_phase_transformer][AT.i_to][:, 0, 0]
    assert np.all(i_from[[1, 3]] == 0.0)
    assert np.all(i_to[[2, 3]] == 0.0)
    assert np.all(i_from[[0, 2, 4]] > 0.0)
    assert np.all(batch_result[CT.node][AT.u][[1, 2, 3], 1] == 0.0)
    assert np.all(batch_result[CT.node][AT.u][[0, 4], 1, 0] > 0.0)


def _split_phase_input():
    """node 0 (abc, 22 kV, source) -- single_phase_transformer 400 a -> ab_center_tapped (12.7 kV / 240-0-240 V,
    solidly earthed centre tap) -- node 1 (split phase, u_rated 240 sqrt3) with unbalanced loads on the legs a and b"""
    input_data = _single_phase_transformer_input()
    input_data[CT.node][AT.u_rated] = [22e3, 240.0 * np.sqrt(3.0)]
    transformer_input = input_data[CT.single_phase_transformer]
    transformer_input[AT.u1] = [22e3 / np.sqrt(3.0)]
    transformer_input[AT.u2] = [480.0]
    transformer_input[AT.sn] = [25e3]
    transformer_input[AT.pk] = [400.0]
    transformer_input[AT.from_connection] = [WindingConnection1P.a]
    transformer_input[AT.to_connection] = [WindingConnection1P.ab_center_tapped]
    transformer_input[AT.r_grounding_to] = [np.nan]
    input_data[CT.asym_load][AT.p_specified] = [[8e3, 2e3, 0.0]]
    input_data[CT.asym_load][AT.q_specified] = [[2e3, 0.5e3, 0.0]]
    return input_data


@pytest.mark.parametrize(
    ("calculation_method", "calculation_initialization"),
    [
        (CalculationMethod.newton_raphson, CalculationInitialization.default),
        (CalculationMethod.newton_raphson, CalculationInitialization.flat),
        (CalculationMethod.newton_raphson, CalculationInitialization.average_source),
        (CalculationMethod.iterative_current, CalculationInitialization.default),
    ],
)
def test_single_phase_transformer_center_tapped(calculation_method, calculation_initialization):
    model = PowerGridModelWithExt(_split_phase_input())
    result = model.calculate_power_flow_with_ext(
        symmetric=False,
        calculation_method=calculation_method,
        calculation_initialization=calculation_initialization,
        experimental_features=_ExperimentalFeatures.enabled,
    )
    u = result[CT.node][AT.u][1]
    u_angle = result[CT.node][AT.u_angle][1]
    # the legs of the split-phase node are close to opposite, leg a is loaded more, and there is no phase c
    assert 0.0 < u[0] < u[1]
    assert u[2] == 0.0
    leg_angle = np.remainder(u_angle[1] - u_angle[0], 2.0 * np.pi)
    assert leg_angle == pytest.approx(np.pi, abs=0.05)

    # each half is rated sn / 2
    transformer = result[CT.single_phase_transformer]
    s_from = np.sum(transformer[AT.p_from][0] + 1j * transformer[AT.q_from][0])
    s_to = transformer[AT.p_to][0] + 1j * transformer[AT.q_to][0]
    assert s_to[2] == 0.0
    expected_loading = max(abs(s_from), 2.0 * abs(s_to[0]), 2.0 * abs(s_to[1])) / 25e3
    assert transformer[AT.loading][0] == pytest.approx(expected_loading, rel=1e-12)
    assert transformer[AT.loading][0] == pytest.approx(2.0 * abs(s_to[0]) / 25e3, rel=1e-12)


def test_single_phase_transformer_center_tapped_invalid_input():
    input_data = _split_phase_input()
    input_data[CT.single_phase_transformer][AT.from_connection] = [WindingConnection1P.ab_center_tapped]
    with pytest.raises(InvalidPhaseConnection, match="only allowed on the to side"):
        PowerGridModel(input_data)

    # uk_half and pk_half only apply to a centre-tapped winding
    for attribute in (AT.uk_half, AT.pk_half):
        input_data = _single_phase_transformer_input()
        input_data[CT.single_phase_transformer][attribute] = [0.03]
        with pytest.raises(InvalidPhaseConnection, match=attribute):
            PowerGridModel(input_data)


def test_single_phase_transformer_center_tapped_short_circuit():
    enabled = _ExperimentalFeatures.enabled
    input_data = _split_phase_input()
    fault_input = initialize_array(DatasetType.input, CT.fault, 1)
    fault_input[AT.id] = [600]
    fault_input[AT.status] = [1]
    fault_input[AT.fault_type] = [FaultType.two_phase]
    fault_input[AT.fault_phase] = [FaultPhase.ab]
    fault_input[AT.fault_object] = [1]
    input_data[CT.fault] = fault_input
    sc_result = PowerGridModelWithExt(input_data).calculate_short_circuit_with_ext(experimental_features=enabled)
    i_f = sc_result[CT.fault][AT.i_f][0]
    assert i_f[0] > 0.0
    assert i_f[1] == pytest.approx(i_f[0])
    assert i_f[2] == 0.0

    # the default phases of a two-phase fault are bc, and the split-phase node has no phase c
    input_data[CT.fault][AT.fault_phase] = [FaultPhase.default_value]
    with pytest.raises(
        InvalidShortCircuitPhaseOrType,
        match=r"of fault 600 are not all present at its node 1 .*default of the fault type",
    ):
        PowerGridModelWithExt(input_data).calculate_short_circuit_with_ext(experimental_features=enabled)


def _fault_phase_input(fault_type: FaultType, fault_phase: FaultPhase, node: int, status: int = 1):
    """node 0 (abc, source) -- line 200 (a) -- node 1 (a); node 0 -- line 201 (b) and line 202 (a) -- node 2 (ab);
    fault 600 of the given type and phases at the given node"""
    input_data = _reduced_phase_input(
        [PhaseConnection.a, PhaseConnection.b, PhaseConnection.a], [(0, 1), (0, 2), (0, 2)], 3
    )
    fault_input = initialize_array(DatasetType.input, CT.fault, 1)
    fault_input[AT.id] = [600]
    fault_input[AT.status] = [status]
    fault_input[AT.fault_type] = [fault_type]
    fault_input[AT.fault_phase] = [fault_phase]
    fault_input[AT.fault_object] = [node]
    input_data[CT.fault] = fault_input
    return input_data


@pytest.mark.parametrize(
    ("fault_type", "fault_phase", "node", "faulted", "valid"),
    [
        pytest.param(FaultType.three_phase, FaultPhase.abc, 0, "abc", True, id="three_phase_at_abc"),
        pytest.param(FaultType.single_phase_to_ground, FaultPhase.c, 0, "c", True, id="spg_c_at_abc"),
        pytest.param(FaultType.two_phase, FaultPhase.default_value, 0, "bc", True, id="two_phase_default_at_abc"),
        pytest.param(FaultType.single_phase_to_ground, FaultPhase.a, 1, "a", True, id="spg_a_at_a"),
        pytest.param(FaultType.single_phase_to_ground, FaultPhase.default_value, 1, "a", True, id="spg_default_at_a"),
        pytest.param(FaultType.single_phase_to_ground, FaultPhase.b, 1, "b", False, id="spg_b_at_a"),
        pytest.param(FaultType.three_phase, FaultPhase.abc, 1, "abc", False, id="three_phase_at_a"),
        pytest.param(FaultType.two_phase, FaultPhase.ab, 1, "ab", False, id="two_phase_ab_at_a"),
        pytest.param(FaultType.two_phase, FaultPhase.default_value, 1, "bc", False, id="two_phase_default_at_a"),
        pytest.param(FaultType.single_phase_to_ground, FaultPhase.b, 2, "b", True, id="spg_b_at_ab"),
        pytest.param(FaultType.two_phase, FaultPhase.ab, 2, "ab", True, id="two_phase_ab_at_ab"),
        pytest.param(FaultType.two_phase_to_ground, FaultPhase.ab, 2, "ab", True, id="two_phase_to_ground_ab_at_ab"),
        pytest.param(FaultType.single_phase_to_ground, FaultPhase.c, 2, "c", False, id="spg_c_at_ab"),
        pytest.param(FaultType.three_phase, FaultPhase.abc, 2, "abc", False, id="three_phase_at_ab"),
        pytest.param(
            FaultType.two_phase_to_ground,
            FaultPhase.default_value,
            2,
            "bc",
            False,
            id="two_phase_to_ground_default_at_ab",
        ),
        pytest.param(FaultType.two_phase, FaultPhase.default_value, 2, "bc", False, id="two_phase_default_at_ab"),
        pytest.param(FaultType.three_phase, FaultPhase.default_value, 1, "abc", False, id="three_phase_default_at_a"),
        pytest.param(FaultType.three_phase, FaultPhase.default_value, 2, "abc", False, id="three_phase_default_at_ab"),
        pytest.param(
            FaultType.two_phase_to_ground,
            FaultPhase.default_value,
            1,
            "bc",
            False,
            id="two_phase_to_ground_default_at_a",
        ),
    ],
)
def test_fault_phases_at_reduced_phase_nodes(fault_type, fault_phase, node, faulted, valid):
    model = PowerGridModelWithExt(_fault_phase_input(fault_type, fault_phase, node))
    enabled = _ExperimentalFeatures.enabled
    if not valid:
        with pytest.raises(InvalidShortCircuitPhaseOrType, match="of fault 600 ") as error:
            model.calculate_short_circuit_with_ext(experimental_features=enabled)
        message = str(error.value)
        assert f"phases ({faulted}) of fault 600 are not all present at its node {node} " in message
        # the hint only appears when the phases are the default and an explicit phase of this fault type fits the node
        n_node_phases = {0: 3, 1: 1, 2: 2}[node]
        min_node_phases = {
            FaultType.three_phase: None,  # needs abc, so no explicit phase helps at a node without all phases
            FaultType.single_phase_to_ground: 1,
            FaultType.two_phase: 2,
            FaultType.two_phase_to_ground: 2,
        }[fault_type]
        explicit_phase_can_fit = min_node_phases is not None and n_node_phases >= min_node_phases
        hint = "the phases are the default of the fault type, so set fault_phase explicitly"
        assert (hint in message) == (fault_phase == FaultPhase.default_value and explicit_phase_can_fit)
        if not explicit_phase_can_fit:
            assert "this short circuit type cannot be applied there\n" in message
        return

    i_f = model.calculate_short_circuit_with_ext(experimental_features=enabled)[CT.fault][AT.i_f][0]
    for phase, letter in enumerate("abc"):
        if letter in faulted:
            assert i_f[phase] > 0.0
        else:
            assert i_f[phase] == 0.0


def test_fault_phases_status_off_and_isolated():
    enabled = _ExperimentalFeatures.enabled
    # a fault with status off is not part of the calculation
    model = PowerGridModelWithExt(_fault_phase_input(FaultType.single_phase_to_ground, FaultPhase.b, 1, status=0))
    i_f = model.calculate_short_circuit_with_ext(experimental_features=enabled)[CT.fault][AT.i_f][0]
    assert i_f.tolist() == [0.0, 0.0, 0.0]

    # a fault at an isolated node has no current
    input_data = _fault_phase_input(FaultType.single_phase_to_ground, FaultPhase.b, 1)
    input_data[CT.asym_line][AT.from_status][0] = 0
    model = PowerGridModelWithExt(input_data)
    i_f = model.calculate_short_circuit_with_ext(experimental_features=enabled)[CT.fault][AT.i_f][0]
    assert i_f.tolist() == [0.0, 0.0, 0.0]


def test_fault_phases_batch_fixed_in_every_scenario():
    # the cache run with the base input fails, but every scenario provides a valid fault
    enabled = _ExperimentalFeatures.enabled
    model = PowerGridModelWithExt(_fault_phase_input(FaultType.single_phase_to_ground, FaultPhase.b, 1))

    update = initialize_array(DatasetType.update, CT.fault, (2, 1))
    update[AT.id] = 600
    update[AT.fault_phase] = [[FaultPhase.a], [FaultPhase.b]]
    update[AT.fault_object] = [[1], [2]]
    result = model.calculate_short_circuit_with_ext(update_data={CT.fault: update}, experimental_features=enabled)
    i_f = result[CT.fault][AT.i_f]
    assert i_f[0, 0, 0] > 0.0
    assert i_f[1, 0, 1] > 0.0

    # a single, non-batch calculation on the same model still raises
    with pytest.raises(InvalidShortCircuitPhaseOrType, match="of fault 600 are not all present at its node 1 "):
        model.calculate_short_circuit_with_ext(experimental_features=enabled)


def test_fault_phases_batch_scenario_leaves_invalid_base_fault():
    # the base fault b at the a node 1 is invalid; scenario 0 fixes it, scenario 1 leaves it unchanged and fails on its
    # own with the fault ID
    enabled = _ExperimentalFeatures.enabled
    model = PowerGridModelWithExt(_fault_phase_input(FaultType.single_phase_to_ground, FaultPhase.b, 1))

    update = initialize_array(DatasetType.update, CT.fault, (2, 1))
    update[AT.id] = 600
    update[AT.fault_phase][0] = FaultPhase.a  # scenario 1 keeps the null value: nothing to update
    result = model.calculate_short_circuit_with_ext(
        update_data={CT.fault: update}, continue_on_batch_error=True, experimental_features=enabled
    )
    assert result[CT.fault][AT.i_f][0, 0, 0] > 0.0
    assert model.batch_error is not None
    assert model.batch_error.failed_scenarios.tolist() == [1]
    assert model.batch_error.errors[0].__class__ is InvalidShortCircuitPhaseOrType
    assert "phases (b) of fault 600 are not all present at its node 1 " in model.batch_error.error_messages[0]

    # without continue_on_batch_error, the batch raises
    with pytest.raises(PowerGridBatchError):
        model.calculate_short_circuit_with_ext(update_data={CT.fault: update}, experimental_features=enabled)


def test_fault_phases_on_node_merged_by_link():
    # node 3 is connected to the a node 1 by a link: both are one bus with phase a, and the error names the user node
    enabled = _ExperimentalFeatures.enabled
    input_data = _reduced_phase_input([PhaseConnection.a], [(0, 1)], 4)
    link_input = initialize_array(DatasetType.input, CT.link, 1)
    link_input[AT.id] = [400]
    link_input[AT.from_node] = [1]
    link_input[AT.to_node] = [3]
    link_input[AT.from_status] = [1]
    link_input[AT.to_status] = [1]
    input_data[CT.link] = link_input
    fault_input = initialize_array(DatasetType.input, CT.fault, 1)
    fault_input[AT.id] = [600]
    fault_input[AT.status] = [1]
    fault_input[AT.fault_type] = [FaultType.single_phase_to_ground]
    fault_input[AT.fault_phase] = [FaultPhase.a]
    fault_input[AT.fault_object] = [3]
    input_data[CT.fault] = fault_input
    model = PowerGridModelWithExt(input_data)

    i_f = model.calculate_short_circuit_with_ext(experimental_features=enabled)[CT.fault][AT.i_f][0]
    assert i_f[0] > 0.0
    assert i_f[1] == 0.0
    assert i_f[2] == 0.0

    model.update(update_data={CT.fault: _fault_update(FaultPhase.b)})
    with pytest.raises(InvalidShortCircuitPhaseOrType) as error:
        model.calculate_short_circuit_with_ext(experimental_features=enabled)
    assert "phases (b) of fault 600 are not all present at its node 3 (node phases: a)" in str(error.value)


def _fault_update(fault_phase: FaultPhase):
    update = initialize_array(DatasetType.update, CT.fault, 1)
    update[AT.id] = 600
    update[AT.fault_phase] = fault_phase
    return update


def test_fault_phases_batch_topology_change():
    # the base fault b at the ab node 2 is valid; switching off the b line makes node 2 an a node in scenario 1
    enabled = _ExperimentalFeatures.enabled
    model = PowerGridModelWithExt(_fault_phase_input(FaultType.single_phase_to_ground, FaultPhase.b, 2))

    update = initialize_array(DatasetType.update, CT.asym_line, (2, 1))
    update[AT.id] = 201
    update[AT.from_status] = [[1], [0]]
    update[AT.to_status] = [[1], [0]]
    result = model.calculate_short_circuit_with_ext(
        update_data={CT.asym_line: update}, continue_on_batch_error=True, experimental_features=enabled
    )
    assert result[CT.fault][AT.i_f][0, 0, 1] > 0.0
    assert model.batch_error is not None
    assert model.batch_error.failed_scenarios.tolist() == [1]
    assert model.batch_error.errors[0].__class__ is InvalidShortCircuitPhaseOrType
    assert "fault 600" in model.batch_error.error_messages[0]
    assert "(node phases: a)" in model.batch_error.error_messages[0]


def test_load_gen_connection_exported():
    assert [(member.name, member.value) for member in LoadGenConnection] == [("wye", 0), ("delta", 1)]


def _delta_input(
    p,
    q,
    *,
    load_type: CT = CT.asym_load,
    connection: int | None = LoadGenConnection.delta,
    load_gen_type: LoadGenType = LoadGenType.const_power,
    line_phases: PhaseConnection | None = None,
):
    """
    node 0 (abc, 400 V, stiff source) with a load of the given type and connection at node 0, or, if line_phases is
    given, behind an asym_line 200 with those phases at node 1 (which then has the phases of the line)
    """
    n_node = 1 if line_phases is None else 2
    node_input = initialize_array(DatasetType.input, CT.node, n_node)
    node_input[AT.id] = list(range(n_node))
    node_input[AT.u_rated] = 400.0

    source_input = initialize_array(DatasetType.input, CT.source, 1)
    source_input[AT.id] = [100]
    source_input[AT.node] = [0]
    source_input[AT.status] = [1]
    source_input[AT.u_ref] = [1.0]
    source_input[AT.sk] = [1e12]

    load_input = initialize_array(DatasetType.input, load_type, 1)
    load_input[AT.id] = [300]
    load_input[AT.node] = [n_node - 1]
    load_input[AT.status] = [1]
    load_input[AT.type] = [load_gen_type]
    load_input[AT.p_specified] = [p]
    load_input[AT.q_specified] = [q]
    if connection is not None:
        load_input[AT.connection] = [connection]

    input_data = {CT.node: node_input, CT.source: source_input, load_type: load_input}

    if line_phases is not None:
        line_input = initialize_array(DatasetType.input, CT.asym_line, 1)
        line_input[AT.id] = [200]
        line_input[AT.from_node] = [0]
        line_input[AT.to_node] = [1]
        line_input[AT.from_status] = [1]
        line_input[AT.to_status] = [1]
        line_input["phases"] = [line_phases]
        for name in ("aa",) if line_phases == PhaseConnection.a else ("aa", "bb", "ba"):
            line_input[f"r_{name}"] = [0.5 if name[0] == name[1] else 0.1]
            line_input[f"x_{name}"] = [0.4 if name[0] == name[1] else 0.1]
            line_input[f"c_{name}"] = [1e-8 if name[0] == name[1] else 0.0]
        input_data[CT.asym_line] = line_input

    return input_data


@pytest.mark.parametrize("method", [CalculationMethod.newton_raphson, CalculationMethod.iterative_current])
@pytest.mark.parametrize(
    "load_gen_type", [LoadGenType.const_power, LoadGenType.const_current, LoadGenType.const_impedance]
)
def test_delta_load_closed_form(method, load_gen_type):
    # a stiff source: the ab branch of the load sees 400 V, so p/q_specified is the power of the branch
    p = 1e4
    q = 0.0
    model = PowerGridModelWithExt(_delta_input([p, 0.0, 0.0], [q, 0.0, 0.0], load_gen_type=load_gen_type))
    result = model.calculate_power_flow_with_ext(
        symmetric=False, calculation_method=method, experimental_features=_ExperimentalFeatures.enabled
    )[CT.asym_load]

    # a resistive ab load is shared over a and b, with opposite reactive power; c is not touched
    assert result[AT.p][0] == pytest.approx([p / 2, p / 2, 0.0], rel=1e-5, abs=1e-6)
    assert result[AT.q][0] == pytest.approx([-p / (2 * np.sqrt(3.0)), p / (2 * np.sqrt(3.0)), 0.0], rel=1e-5, abs=1e-6)
    assert result[AT.i][0] == pytest.approx([p / 400.0, p / 400.0, 0.0], rel=1e-5, abs=1e-6)
    assert result[AT.p][0][2] == 0.0
    assert result[AT.q][0][2] == 0.0
    assert result[AT.i][0][2] == 0.0


def test_delta_gen_untouched_phase_is_positive_zero():
    # an asym_load reports -0 for a zero power, because its output is the negated injection; an asym_gen does not
    model = PowerGridModelWithExt(_delta_input([-1e4, 0.0, 0.0], [-1e3, 0.0, 0.0], load_type=CT.asym_gen))
    result = model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)
    assert result[CT.asym_gen][AT.p][0].sum() == pytest.approx(-1e4, rel=1e-5)
    assert result[CT.asym_gen][AT.p][0][2] == 0.0
    assert not np.signbit(result[CT.asym_gen][AT.p][0][2])
    assert not np.signbit(result[CT.asym_gen][AT.q][0][2])
    assert not np.signbit(result[CT.asym_gen][AT.i][0][2])


def test_delta_load_total_power():
    model = PowerGridModelWithExt(_delta_input([3e3, -1e3, 2e3], [1e3, 5e2, -2e3]))
    result = model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)
    # the total power of the three branches is the sum of the phase powers
    assert result[CT.asym_load][AT.p][0].sum() == pytest.approx(4e3, rel=1e-6)
    assert result[CT.asym_load][AT.q][0].sum() == pytest.approx(-5e2, rel=1e-5)
    assert result[CT.source][AT.p][0].sum() == pytest.approx(4e3, rel=1e-3)


def test_delta_load_experimental_gate():
    model = PowerGridModelWithExt(_delta_input([1e4, 0.0, 0.0], [0.0, 0.0, 0.0]))
    with pytest.raises(ExperimentalFeature, match="delta-connected load or generator"):
        model.calculate_power_flow(symmetric=False)
    with pytest.raises(ExperimentalFeature):
        model.calculate_short_circuit()


def test_delta_load_unsupported_calculations():
    model = PowerGridModelWithExt(_delta_input([1e4, 0.0, 0.0], [0.0, 0.0, 0.0]))
    enabled = _ExperimentalFeatures.enabled

    with pytest.raises(ReducedPhaseNotSupported, match="delta-connected loads or generators"):
        model.calculate_power_flow_with_ext(symmetric=True, experimental_features=enabled)
    for symmetric in (True, False):
        with pytest.raises(ReducedPhaseNotSupported, match="delta-connected loads or generators"):
            model.calculate_state_estimation_with_ext(symmetric=symmetric, experimental_features=enabled)


def test_delta_load_short_circuit():
    # the load is ignored in a short circuit calculation
    model = PowerGridModelWithExt(_delta_input([1e4, 0.0, 0.0], [0.0, 0.0, 0.0]))
    result = model.calculate_short_circuit_with_ext(experimental_features=_ExperimentalFeatures.enabled)
    assert result[CT.node][AT.u_pu].shape == (1, 3)


@pytest.mark.parametrize("connection", [None, LoadGenConnection.wye])
def test_wye_connection_is_the_default(connection):
    p = [1e3, 2e3, 3e3]
    q = [1e2, 2e2, 3e2]
    model = PowerGridModelWithExt(_delta_input(p, q, connection=connection))
    # a wye load needs no experimental features and has the specified power per phase
    result = model.calculate_power_flow(symmetric=False)
    assert result[CT.asym_load][AT.p][0] == pytest.approx(p, rel=1e-6)
    assert result[CT.asym_load][AT.q][0] == pytest.approx(q, rel=1e-6)


def test_delta_sym_load_not_supported():
    with pytest.raises(InvalidPhaseConnection, match="delta"):
        PowerGridModel(_delta_input(1e4, 0.0, load_type=CT.sym_load))
    # wye and not provided are fine
    PowerGridModel(_delta_input(1e4, 0.0, load_type=CT.sym_load, connection=LoadGenConnection.wye))
    PowerGridModel(_delta_input(1e4, 0.0, load_type=CT.sym_load, connection=None))


def test_load_gen_invalid_connection():
    with pytest.raises(MissingCaseForEnumError):
        PowerGridModel(_delta_input([1e4, 0.0, 0.0], [0.0, 0.0, 0.0], connection=2))


@pytest.mark.parametrize(
    ("line_phases", "p", "valid"),
    [
        pytest.param(PhaseConnection.a, [1e3, 0.0, 0.0], False, id="a_node_ab"),
        pytest.param(PhaseConnection.a, [0.0, 0.0, 1e3], False, id="a_node_ca"),
        pytest.param(PhaseConnection.ab, [1e3, 0.0, 0.0], True, id="ab_node_ab"),
        pytest.param(PhaseConnection.ab, [0.0, 1e3, 0.0], False, id="ab_node_bc"),
        pytest.param(PhaseConnection.ab, [0.0, 0.0, 1e3], False, id="ab_node_ca"),
        pytest.param(PhaseConnection.ab, [1e3, 0.0, 0.0], True, id="ab_node_ab_other_entries_zero"),
        pytest.param(PhaseConnection.ab, [1e3, np.nan, np.nan], True, id="ab_node_other_entries_nan"),
    ],
)
def test_delta_load_on_absent_phase(line_phases, p, valid):
    model = PowerGridModelWithExt(_delta_input(p, [0.0, 0.0, 0.0], line_phases=line_phases))
    enabled = _ExperimentalFeatures.enabled
    if valid:
        result = model.calculate_power_flow_with_ext(symmetric=False, experimental_features=enabled)
        assert result[CT.asym_load][AT.p][0][0] > 0.0
        assert result[CT.asym_load][AT.p][0][2] == 0.0
    else:
        with pytest.raises(LoadGenOnAbsentPhase, match=r"Load or generator 300 .*A delta \(phase-to-phase\) entry"):
            model.calculate_power_flow_with_ext(symmetric=False, experimental_features=enabled)


def test_delta_load_batch_update_matches_separate_runs():
    enabled = _ExperimentalFeatures.enabled
    model = PowerGridModelWithExt(_delta_input([1e4, 0.0, 0.0], [0.0, 0.0, 0.0], line_phases=PhaseConnection.ab))

    p_scenarios = [[1e4, 0.0, 0.0], [5e3, 0.0, 0.0], [2e3, 0.0, 0.0]]
    update = initialize_array(DatasetType.update, CT.asym_load, (3, 1))
    update[AT.id] = 300
    update[AT.p_specified] = [[p] for p in p_scenarios]
    batch = model.calculate_power_flow_with_ext(
        symmetric=False, update_data={CT.asym_load: update}, experimental_features=enabled
    )

    for k, p in enumerate(p_scenarios):
        single = PowerGridModelWithExt(_delta_input(p, [0.0, 0.0, 0.0], line_phases=PhaseConnection.ab))
        expected = single.calculate_power_flow_with_ext(symmetric=False, experimental_features=enabled)
        for component, attribute in ((CT.asym_load, AT.p), (CT.asym_load, AT.q), (CT.node, AT.u)):
            assert batch[component][attribute][k] == pytest.approx(expected[component][attribute], rel=1e-8, abs=1e-8)
        assert batch[CT.asym_load][AT.p][k, 0, 0] > 0.0


def test_delta_gen_with_voltage_regulator_not_supported():
    input_data = _delta_input([-1e4, 0.0, 0.0], [0.0, 0.0, 0.0], load_type=CT.asym_gen)
    regulator_input = initialize_array(DatasetType.input, CT.voltage_regulator, 1)
    regulator_input[AT.id] = [400]
    regulator_input[AT.regulated_object] = [300]
    regulator_input[AT.status] = [1]
    regulator_input[AT.u_ref] = [1.0]
    input_data[CT.voltage_regulator] = regulator_input
    model = PowerGridModelWithExt(input_data)
    with pytest.raises(UnsupportedLoadGenTypeForVoltageRegulator):
        model.calculate_power_flow_with_ext(symmetric=False, experimental_features=_ExperimentalFeatures.enabled)
