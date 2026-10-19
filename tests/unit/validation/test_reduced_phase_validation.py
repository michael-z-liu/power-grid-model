# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""
Validation of reduced-phase (one- or two-phase) networks.

The "mirror" tests check that an input flagged by the Python validation is also rejected by the C++ core, so that the
two do not drift apart.
"""

import warnings
from typing import Any

import numpy as np
import pytest

from power_grid_model import (
    AttributeType as AT,
    CalculationMethod,
    ComponentType as CT,
    DatasetType,
    PowerGridModel,
    initialize_array,
)
from power_grid_model.enum import (
    AngleMeasurementType,
    Branch3Side,
    BranchSide,
    CalculationType,
    LoadGenConnection,
    LoadGenType,
    MeasuredTerminalType,
    PhaseConnection,
    TapChangingStrategy,
    WindingConnection1P,
    WindingType,
)
from power_grid_model.errors import (
    AutomaticTapInputError,
    InvalidPhaseConnection,
    InvalidTransformerClock,
    IterationDiverge,
    MissingCaseForEnumError,
    UnsupportedLoadGenTypeForVoltageRegulator,
)
from power_grid_model.validation import (
    ValidationWarning,
    assert_valid_batch_data,
    assert_valid_input_data,
    validate_batch_data,
    validate_input_data,
)
from power_grid_model.validation._validation import _connected_components
from power_grid_model.validation.errors import (
    InvalidAssociatedEnumValueError,
    InvalidEnumValueError,
    InvalidIdError,
    InvalidValueError,
    InvalidVoltageRegulationError,
    MissingValueError,
    MultiFieldValidationError,
    NotApplicableValueError,
    NotBetweenError,
    NotBetweenOrAtError,
    NotGreaterOrEqualError,
    NotGreaterThanError,
    NotLessThanError,
    SinglePhaseTransformerClockError,
)

LINE_ID = 3

# (field, row conductor, column conductor) of the r, x and c matrices; "n" is the neutral
ENTRIES = [
    *[(AT(f"r_{row}{col}"), row, col) for row, col in ("aa", "ba", "bb", "ca", "cb", "cc", "na", "nb", "nc", "nn")],
    *[(AT(f"x_{row}{col}"), row, col) for row, col in ("aa", "ba", "bb", "ca", "cb", "cc", "na", "nb", "nc", "nn")],
    *[(AT(f"c_{row}{col}"), row, col) for row, col in ("aa", "ba", "bb", "ca", "cb", "cc")],
]
SELF_R, MUTUAL_R, SELF_R_NEUTRAL, MUTUAL_R_NEUTRAL = 1.0, 0.1, 0.5, 0.05
SELF_C, MUTUAL_C = 1e-8, -1e-9
REDUCED_PHASES = {
    PhaseConnection.a: "a",
    PhaseConnection.b: "b",
    PhaseConnection.c: "c",
    PhaseConnection.ab: "ab",
    PhaseConnection.ac: "ac",
    PhaseConnection.bc: "bc",
}


def _entry_value(field: AT, row: str, col: str) -> float:
    if field.startswith("c_"):
        return SELF_C if row == col else MUTUAL_C
    if row == col:
        return SELF_R_NEUTRAL if row == "n" else SELF_R
    return MUTUAL_R_NEUTRAL if row == "n" else MUTUAL_R


def _present(phases: PhaseConnection, with_neutral: bool) -> set[str]:
    return set(REDUCED_PHASES[phases]) | ({"n"} if with_neutral else set())


def _swer_asym_line(phases: PhaseConnection, with_neutral: bool = False) -> np.ndarray:
    """A one- or two-phase asym_line with exactly the entries that its phases (and neutral) need."""
    asym_line = initialize_array(DatasetType.input, CT.asym_line, 1)
    asym_line[AT.id] = LINE_ID
    asym_line[AT.from_node] = 1
    asym_line[AT.to_node] = 2
    asym_line[AT.from_status] = 1
    asym_line[AT.to_status] = 1
    asym_line[AT.phases] = phases
    present = _present(phases, with_neutral)
    for field, row, col in ENTRIES:
        if row in present and col in present:
            asym_line[field] = _entry_value(field, row, col)
    return asym_line


def _input(asym_line: np.ndarray) -> dict[CT, np.ndarray]:
    node = initialize_array(DatasetType.input, CT.node, 2)
    node[AT.id] = [1, 2]
    node[AT.u_rated] = 22e3
    return {CT.node: node, CT.asym_line: asym_line}


def _validate(asym_line: np.ndarray):
    return validate_input_data(_input(asym_line), symmetric=False)


def _assert_cpp_rejects(asym_line: np.ndarray):
    with pytest.raises(InvalidPhaseConnection):
        PowerGridModel(_input(asym_line))


def _configurations():
    for phases in REDUCED_PHASES:
        for with_neutral in (False, True):
            yield phases, with_neutral


def _config_id(phases: PhaseConnection, with_neutral: bool) -> str:
    return f"{phases.name}{'-n' if with_neutral else ''}"


def _required_cases():
    for phases, with_neutral in _configurations():
        present = _present(phases, with_neutral)
        for field, row, col in ENTRIES:
            if row in present and col in present:
                yield pytest.param(phases, with_neutral, field, id=f"{_config_id(phases, with_neutral)}-{field}")


def _not_applicable_cases():
    for phases, with_neutral in _configurations():
        present = _present(phases, with_neutral)
        for field, row, col in ENTRIES:
            # r_nn or x_nn on a line without a neutral adds the neutral instead
            if not (row in present and col in present) and field not in (AT.r_nn, AT.x_nn):
                yield pytest.param(
                    phases, with_neutral, field, row, col, id=f"{_config_id(phases, with_neutral)}-{field}"
                )


@pytest.mark.parametrize(
    ("phases", "with_neutral"),
    [pytest.param(*config, id=_config_id(*config)) for config in _configurations()],
)
def test_reduced_phase_asym_line_valid(phases: PhaseConnection, with_neutral: bool):
    asym_line = _swer_asym_line(phases, with_neutral)
    assert _validate(asym_line) is None
    PowerGridModel(_input(asym_line))


@pytest.mark.parametrize(("phases", "with_neutral", "field"), list(_required_cases()))
def test_reduced_phase_asym_line_missing_entry(phases: PhaseConnection, with_neutral: bool, field: AT):
    asym_line = _swer_asym_line(phases, with_neutral)
    asym_line[field] = np.nan
    assert _validate(asym_line) == [MissingValueError(CT.asym_line, field, [LINE_ID])]
    _assert_cpp_rejects(asym_line)


@pytest.mark.parametrize(("phases", "with_neutral", "field", "row", "col"), list(_not_applicable_cases()))
def test_reduced_phase_asym_line_not_applicable_entry(
    phases: PhaseConnection, with_neutral: bool, field: AT, row: str, col: str
):
    asym_line = _swer_asym_line(phases, with_neutral)
    asym_line[field] = _entry_value(field, row, col)
    assert _validate(asym_line) == [NotApplicableValueError(CT.asym_line, [AT.phases, field], [LINE_ID])]
    _assert_cpp_rejects(asym_line)


def test_reduced_phase_asym_line_cases_cover_all_entries():
    required = {field for _, _, field in (case.values for case in _required_cases())}
    not_applicable = {field for _, _, field, _, _ in (case.values for case in _not_applicable_cases())}
    all_fields = {field for field, _, _ in ENTRIES}
    assert required == all_fields
    # r_nn and x_nn define the neutral, so they always apply
    assert not_applicable == all_fields - {AT.r_nn, AT.x_nn}


def test_two_phase_asym_line_neutral():
    # ab line with a neutral: r/x_nn and r/x_na, r/x_nb
    asym_line = _swer_asym_line(PhaseConnection.ab, with_neutral=True)
    assert asym_line[AT.r_nn] == SELF_R_NEUTRAL
    assert asym_line[AT.r_na] == MUTUAL_R_NEUTRAL
    assert asym_line[AT.r_nb] == MUTUAL_R_NEUTRAL
    assert np.isnan(asym_line[AT.r_nc])
    assert _validate(asym_line) is None

    # x_nn missing while r_nn is provided
    missing_x_nn = asym_line.copy()
    missing_x_nn[AT.x_nn] = np.nan
    errors = _validate(missing_x_nn)
    assert errors == [MissingValueError(CT.asym_line, AT.x_nn, [LINE_ID])]
    _assert_cpp_rejects(missing_x_nn)

    # r_nc refers to the absent phase c
    with_r_nc = asym_line.copy()
    with_r_nc[AT.r_nc] = 0.05
    errors = _validate(with_r_nc)
    assert errors == [NotApplicableValueError(CT.asym_line, [AT.phases, AT.r_nc], [LINE_ID])]
    _assert_cpp_rejects(with_r_nc)

    # r_na without a neutral (no r_nn or x_nn)
    no_neutral = _swer_asym_line(PhaseConnection.ab)
    no_neutral[AT.r_na] = 0.05
    errors = _validate(no_neutral)
    assert errors == [NotApplicableValueError(CT.asym_line, [AT.phases, AT.r_na], [LINE_ID])]
    _assert_cpp_rejects(no_neutral)

    # the abc rule "all neutral entries or none" does not apply
    for errors in (_validate(missing_x_nn), _validate(with_r_nc), _validate(no_neutral)):
        assert not [error for error in errors or [] if type(error) is MultiFieldValidationError]


def test_one_phase_asym_line_c0_c1():
    asym_line = _swer_asym_line(PhaseConnection.a)
    asym_line[AT.c0] = 1e-8
    asym_line[AT.c1] = 1e-8
    errors = _validate(asym_line)
    assert errors is not None
    assert sorted(errors, key=str) == [
        NotApplicableValueError(CT.asym_line, [AT.phases, AT.c0], [LINE_ID]),
        NotApplicableValueError(CT.asym_line, [AT.phases, AT.c1], [LINE_ID]),
    ]
    assert not [error for error in errors if type(error) is MultiFieldValidationError]
    _assert_cpp_rejects(asym_line)


def test_asym_line_unknown_phases():
    asym_line = _swer_asym_line(PhaseConnection.a)
    asym_line[AT.phases] = 7
    for field, _, _ in ENTRIES:
        asym_line[field] = np.nan
    errors = _validate(asym_line)
    assert errors == [InvalidEnumValueError(CT.asym_line, AT.phases, [LINE_ID], PhaseConnection)]
    with pytest.raises(MissingCaseForEnumError):
        PowerGridModel(_input(asym_line))


def test_asym_line_mixed_abc_and_reduced_rows():
    abc_id, a_id = 4, 5
    asym_line = initialize_array(DatasetType.input, CT.asym_line, 2)
    asym_line[AT.id] = [abc_id, a_id]
    asym_line[AT.from_node] = 1
    asym_line[AT.to_node] = 2
    asym_line[AT.from_status] = 1
    asym_line[AT.to_status] = 1
    asym_line[AT.phases] = [PhaseConnection.abc, PhaseConnection.a]
    for field, row, col in ENTRIES[:6] + ENTRIES[10:16]:  # the phase r and x entries
        asym_line[field][0] = _entry_value(field, row, col)
    asym_line[AT.c0][0] = 1e-8
    asym_line[AT.c1][0] = 2e-8
    asym_line[AT.r_bb][0] = np.nan
    asym_line[AT.x_aa][1] = 1.0
    asym_line[AT.c_aa][1] = 1e-8

    errors = _validate(asym_line)
    assert errors is not None
    assert sorted(errors, key=str) == [
        MissingValueError(CT.asym_line, AT.r_aa, [a_id]),
        MissingValueError(CT.asym_line, AT.r_bb, [abc_id]),
    ]


def test_one_phase_asym_line_zero_self_capacitance():
    asym_line = _swer_asym_line(PhaseConnection.a)
    asym_line[AT.c_aa] = 0.0
    assert _validate(asym_line) == [NotGreaterThanError(CT.asym_line, AT.c_aa, [LINE_ID], 0)]


def test_reduced_phase_asym_line_batch():
    asym_line = _swer_asym_line(PhaseConnection.a)
    update = initialize_array(DatasetType.update, CT.asym_line, (2, 1))
    update[AT.id] = LINE_ID
    update[AT.to_status] = [[0], [1]]
    assert validate_batch_data(_input(asym_line), {CT.asym_line: update}, symmetric=False) is None


# single_phase_transformer

SPT_ID = 20
FROM_NODE, TO_NODE = 11, 12
U_MV = 22e3
U_SWER = 12.7e3
U_SPLIT_PHASE = 240.0 * np.sqrt(3.0)
SN, UK, PK = 50e3, 0.035, 500.0
SPT_REQUIRED = {
    AT.id,
    AT.from_node,
    AT.to_node,
    AT.from_status,
    AT.to_status,
    AT.u1,
    AT.u2,
    AT.sn,
    AT.uk,
    AT.pk,
    AT.i0,
    AT.p0,
    AT.from_connection,
    AT.to_connection,
}


def _spt(**overrides) -> np.ndarray:
    """A single_phase_transformer from an ab winding (22 kV) to an a winding (12.7 kV), with only the required fields"""
    spt = initialize_array(DatasetType.input, CT.single_phase_transformer, 1)
    values: dict[str, float] = {
        AT.id: SPT_ID,
        AT.from_node: FROM_NODE,
        AT.to_node: TO_NODE,
        AT.from_status: 1,
        AT.to_status: 1,
        AT.u1: U_MV,
        AT.u2: U_SWER,
        AT.sn: SN,
        AT.uk: UK,
        AT.pk: PK,
        AT.i0: 0.01,
        AT.p0: 100.0,
        AT.from_connection: WindingConnection1P.ab,
        AT.to_connection: WindingConnection1P.a,
    }
    values.update(overrides)
    for field, value in values.items():
        spt[field] = value
    return spt


def _center_tapped_spt(**overrides) -> np.ndarray:
    """A split-phase transformer: 240-0-240 V, so u2 = 480 V"""
    return _spt(**{AT.u2: 480.0, AT.to_connection: WindingConnection1P.ab_center_tapped, **overrides})


def _spt_input(
    spt: np.ndarray, u_rated_to: float | None = None, components: dict[CT, np.ndarray] | None = None
) -> dict[CT, np.ndarray]:
    if u_rated_to is None:
        center_tapped = spt[AT.to_connection][0] == WindingConnection1P.ab_center_tapped
        u_rated_to = U_SPLIT_PHASE if center_tapped else U_MV
    node = initialize_array(DatasetType.input, CT.node, 2)
    node[AT.id] = [FROM_NODE, TO_NODE]
    node[AT.u_rated] = [U_MV, u_rated_to]
    return {CT.node: node, CT.single_phase_transformer: spt, **(components or {})}


def _validate_spt(spt: np.ndarray, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("error", ValidationWarning)
        return validate_input_data(_spt_input(spt, **kwargs), symmetric=False)


def _assert_cpp_raises(spt: np.ndarray, error: type[Exception] = InvalidPhaseConnection):
    with pytest.raises(error):
        PowerGridModel(_spt_input(spt))


def test_single_phase_transformer_required_fields():
    spt = initialize_array(DatasetType.input, CT.single_phase_transformer, 1)
    errors = validate_input_data(_spt_input(spt, u_rated_to=U_MV), symmetric=False)
    assert errors is not None
    missing = [error for error in errors if isinstance(error, MissingValueError)]
    assert all(error.component == CT.single_phase_transformer for error in missing)
    # one error per field; the clock, the tap fields, uk_half/pk_half and the earthing impedances are optional
    assert len(missing) == len(SPT_REQUIRED)
    missing_fields = {error.field for error in missing}
    assert missing_fields == SPT_REQUIRED


def test_single_phase_transformer_defaults_valid():
    spt = _spt()
    assert _validate_spt(spt) is None
    PowerGridModel(_spt_input(spt))


@pytest.mark.parametrize(
    ("field", "value", "expected"),
    [
        (AT.u1, 0.0, NotGreaterThanError(CT.single_phase_transformer, AT.u1, [SPT_ID], 0)),
        (AT.u2, -1.0, NotGreaterThanError(CT.single_phase_transformer, AT.u2, [SPT_ID], 0)),
        (AT.sn, 0.0, NotGreaterThanError(CT.single_phase_transformer, AT.sn, [SPT_ID], 0)),
        (AT.uk, 1.0, NotBetweenError(CT.single_phase_transformer, AT.uk, [SPT_ID], (0, 1))),
        (AT.uk, 0.005, NotGreaterOrEqualError(CT.single_phase_transformer, AT.uk, [SPT_ID], f"{AT.pk}/{AT.sn}")),
        (AT.pk, -1.0, NotGreaterOrEqualError(CT.single_phase_transformer, AT.pk, [SPT_ID], 0)),
        (AT.i0, 1.0, NotLessThanError(CT.single_phase_transformer, AT.i0, [SPT_ID], 1)),
        (AT.i0, 0.001, NotGreaterOrEqualError(CT.single_phase_transformer, AT.i0, [SPT_ID], f"{AT.p0}/{AT.sn}")),
        (AT.p0, -1.0, NotGreaterOrEqualError(CT.single_phase_transformer, AT.p0, [SPT_ID], 0)),
        (AT.uk_min, 1.5, NotBetweenError(CT.single_phase_transformer, AT.uk_min, [SPT_ID], (0, 1))),
        (AT.pk_max, -1.0, NotGreaterOrEqualError(CT.single_phase_transformer, AT.pk_max, [SPT_ID], 0)),
        (AT.tap_size, -1.0, NotGreaterOrEqualError(CT.single_phase_transformer, AT.tap_size, [SPT_ID], 0)),
    ],
)
def test_single_phase_transformer_ranges(field: AT, value: float, expected):
    # Python-only rules, as for transformer
    # the uk/pk limits are only used (and checked) together with their tap limit
    tap_limits = {AT.tap_min: -1, AT.tap_max: 1} if field in (AT.uk_min, AT.pk_max) else {}
    errors = _validate_spt(_spt(**{field: value, **tap_limits}))
    assert errors is not None
    assert expected in errors


@pytest.mark.parametrize(
    ("limit", "tap_limit", "other_tap_limit", "value", "expected"),
    [
        (
            AT.uk_min,
            AT.tap_min,
            AT.tap_max,
            1.5,
            NotBetweenError(CT.single_phase_transformer, AT.uk_min, [SPT_ID], (0, 1)),
        ),
        (
            AT.uk_max,
            AT.tap_max,
            AT.tap_min,
            1.5,
            NotBetweenError(CT.single_phase_transformer, AT.uk_max, [SPT_ID], (0, 1)),
        ),
        (
            AT.pk_min,
            AT.tap_min,
            AT.tap_max,
            -1.0,
            NotGreaterOrEqualError(CT.single_phase_transformer, AT.pk_min, [SPT_ID], 0),
        ),
        (
            AT.pk_max,
            AT.tap_max,
            AT.tap_min,
            -1.0,
            NotGreaterOrEqualError(CT.single_phase_transformer, AT.pk_max, [SPT_ID], 0),
        ),
    ],
)
def test_single_phase_transformer_uk_pk_limit_needs_tap_limit(limit, tap_limit, other_tap_limit, value, expected):
    # without its tap limit, the C++ core ignores the uk/pk limit, so it is not an error
    spt = _spt(**{limit: value})
    assert _validate_spt(spt) is None
    PowerGridModel(_spt_input(spt))
    # a provided tap limit of the other side does not enable the check
    spt = _spt(**{limit: value, other_tap_limit: 1 if other_tap_limit == AT.tap_max else -1})
    assert _validate_spt(spt) is None
    PowerGridModel(_spt_input(spt))
    # with its tap limit, the uk/pk limit is checked
    errors = _validate_spt(_spt(**{limit: value, tap_limit: 1 if tap_limit == AT.tap_max else -1}))
    assert errors is not None
    assert expected in errors


@pytest.mark.parametrize("connection", [0, 8, 9])
@pytest.mark.parametrize("field", [AT.from_connection, AT.to_connection])
def test_single_phase_transformer_unknown_connection(field: AT, connection: int):
    spt = _spt(**{field: connection})
    assert _validate_spt(spt) == [
        InvalidEnumValueError(CT.single_phase_transformer, field, [SPT_ID], WindingConnection1P)
    ]
    _assert_cpp_raises(spt, MissingCaseForEnumError)


def test_single_phase_transformer_center_tapped_from_side():
    spt = _spt(**{AT.from_connection: WindingConnection1P.ab_center_tapped})
    assert _validate_spt(spt) == [
        InvalidValueError(
            CT.single_phase_transformer,
            AT.from_connection,
            [SPT_ID],
            [
                WindingConnection1P.a,
                WindingConnection1P.b,
                WindingConnection1P.c,
                WindingConnection1P.ab,
                WindingConnection1P.ac,
                WindingConnection1P.bc,
            ],
        )
    ]
    _assert_cpp_raises(spt)


def test_single_phase_transformer_center_tapped_to_side():
    spt = _center_tapped_spt()
    assert _validate_spt(spt) is None
    PowerGridModel(_spt_input(spt))


def test_single_phase_transformer_invalid_clock():
    spt = _spt(**{AT.clock: 3})
    assert _validate_spt(spt) == [SinglePhaseTransformerClockError(CT.single_phase_transformer, AT.clock, [SPT_ID])]
    _assert_cpp_raises(spt, InvalidTransformerClock)


@pytest.mark.parametrize("clock", [0, 6, -6, 12, 18])
def test_single_phase_transformer_valid_clock(clock: int):
    spt = _spt(**{AT.clock: clock})
    assert _validate_spt(spt) is None
    PowerGridModel(_spt_input(spt))


def test_single_phase_transformer_grounding_on_phase_to_phase_from_side():
    spt = _spt(**{AT.r_grounding_from: 1.0})
    assert _validate_spt(spt) == [
        NotApplicableValueError(CT.single_phase_transformer, [AT.from_connection, AT.r_grounding_from], [SPT_ID])
    ]
    _assert_cpp_raises(spt)


def test_single_phase_transformer_grounding_on_phase_to_phase_to_side():
    spt = _spt(**{AT.u2: U_MV, AT.to_connection: WindingConnection1P.bc, AT.x_grounding_to: 1.0})
    assert _validate_spt(spt) == [
        NotApplicableValueError(CT.single_phase_transformer, [AT.to_connection, AT.x_grounding_to], [SPT_ID])
    ]
    _assert_cpp_raises(spt)


def test_single_phase_transformer_grounding_valid():
    phase_to_earth = _spt(
        **{AT.from_connection: WindingConnection1P.a, AT.u1: U_SWER, AT.r_grounding_from: 1.0, AT.x_grounding_from: 2.0}
    )
    center_tapped = _center_tapped_spt(**{AT.r_grounding_to: 1.0, AT.x_grounding_to: 2.0})
    for spt in (phase_to_earth, center_tapped):
        assert _validate_spt(spt) is None
        PowerGridModel(_spt_input(spt))


def test_single_phase_transformer_negative_grounding_resistance():
    # Python-only rule
    spt = _spt(**{AT.r_grounding_to: -1.0})
    assert _validate_spt(spt) == [NotGreaterOrEqualError(CT.single_phase_transformer, AT.r_grounding_to, [SPT_ID], 0)]


@pytest.mark.parametrize("field", [AT.uk_half, AT.pk_half])
def test_single_phase_transformer_half_winding_not_applicable(field: AT):
    spt = _spt(**{field: 0.02 if field == AT.uk_half else 200.0})
    assert _validate_spt(spt) == [
        NotApplicableValueError(CT.single_phase_transformer, [AT.to_connection, field], [SPT_ID])
    ]
    _assert_cpp_raises(spt)


def test_single_phase_transformer_half_winding_ranges():
    # Python-only rules (D20)
    assert _validate_spt(_center_tapped_spt(**{AT.uk_half: 0.022, AT.pk_half: 200.0})) is None

    assert _validate_spt(_center_tapped_spt(**{AT.uk_half: UK / 2})) == [
        NotGreaterThanError(CT.single_phase_transformer, AT.uk_half, [SPT_ID], f"{AT.uk}/2")
    ]

    for pk_half in (PK / 4 - 1, PK / 2 + 1):
        assert _validate_spt(_center_tapped_spt(**{AT.pk_half: pk_half})) == [
            NotBetweenOrAtError(CT.single_phase_transformer, AT.pk_half, [SPT_ID], (f"{AT.pk}/4", f"{AT.pk}/2"))
        ]
    for pk_half in (PK / 4, PK / 2):
        assert _validate_spt(_center_tapped_spt(**{AT.pk_half: pk_half})) is None

    # r = 0.03, so pk_half = 750 gives pk_half / (sn / 2) = 0.03 > uk_half = 0.02 > uk / 2
    spt = _center_tapped_spt(**{AT.pk: 1500.0, AT.uk_half: 0.02, AT.pk_half: 750.0})
    assert _validate_spt(spt) == [
        NotGreaterOrEqualError(CT.single_phase_transformer, AT.uk_half, [SPT_ID], f"{AT.pk_half}/({AT.sn}/2)")
    ]


@pytest.mark.parametrize(("uk", "r"), [(0.04, 0.0), (0.04, 0.016), (0.1, 0.1), (0.07, 0.07)])
def test_single_phase_transformer_half_winding_defaults_valid(uk: float, r: float):
    spt = _center_tapped_spt(**{AT.uk: uk, AT.pk: r * SN})
    assert _validate_spt(spt) is None
    PowerGridModel(_spt_input(spt))


def test_single_phase_transformer_half_winding_default_uk_half_too_small():
    # with uk = pk/sn (no reactance), the default uk_half = 0.75 pk/sn is below pk_half/(sn/2) = pk/sn when
    # pk_half = pk/2
    uk = 0.04
    spt = _center_tapped_spt(**{AT.uk: uk, AT.pk: uk * SN, AT.pk_half: uk * SN / 2})
    assert _validate_spt(spt) == [
        NotGreaterOrEqualError(CT.single_phase_transformer, AT.uk_half, [SPT_ID], f"{AT.pk_half}/({AT.sn}/2)")
    ]


def _tap_error(field: AT = AT.tap_pos) -> list[NotBetweenOrAtError]:
    return [NotBetweenOrAtError(CT.single_phase_transformer, field, [SPT_ID], (AT.tap_min, AT.tap_max))]


@pytest.mark.parametrize(
    ("taps", "expected"),
    [
        pytest.param({}, None, id="no-tap-fields"),
        pytest.param({AT.tap_min: -2, AT.tap_max: 2, AT.tap_pos: 2}, None, id="both-at-max"),
        pytest.param({AT.tap_min: -2, AT.tap_max: 2, AT.tap_pos: 3}, _tap_error(), id="both-above-max"),
        pytest.param({AT.tap_min: -2, AT.tap_max: 2, AT.tap_pos: -3}, _tap_error(), id="both-below-min"),
        pytest.param({AT.tap_min: 5, AT.tap_max: -5, AT.tap_pos: 7}, _tap_error(), id="reversed-outside"),
        pytest.param({AT.tap_min: 5, AT.tap_max: -5, AT.tap_pos: 0}, None, id="reversed-inside"),
        pytest.param({AT.tap_min: -2, AT.tap_pos: -3}, _tap_error(), id="only-min-below"),
        pytest.param({AT.tap_min: -2, AT.tap_pos: 50}, None, id="only-min-unbounded-above"),
        pytest.param({AT.tap_nom: 5, AT.tap_min: 9, AT.tap_pos: 10}, _tap_error(), id="only-min-reversed-above"),
        pytest.param({AT.tap_nom: 5, AT.tap_min: 9, AT.tap_pos: 3}, None, id="only-min-reversed-unbounded-below"),
        pytest.param({AT.tap_max: 2, AT.tap_pos: 3}, _tap_error(), id="only-max-above"),
        pytest.param({AT.tap_max: 2, AT.tap_pos: -50}, None, id="only-max-unbounded-below"),
        pytest.param({AT.tap_min: -2, AT.tap_max: 2, AT.tap_nom: 3}, _tap_error() + _tap_error(AT.tap_nom), id="nom"),
        pytest.param(
            {AT.tap_min: -2, AT.tap_max: 2, AT.tap_nom: 3, AT.tap_pos: 0}, _tap_error(AT.tap_nom), id="only-nom"
        ),
    ],
)
def test_single_phase_transformer_taps(taps: dict, expected):
    # Python-only rule; the C++ core clamps tap_pos to the provided limits instead (D18)
    spt = _spt(**taps)
    assert _validate_spt(spt) == expected
    PowerGridModel(_spt_input(spt))


def test_single_phase_transformer_tap_side():
    spt = _spt(**{AT.tap_side: 5})
    assert _validate_spt(spt) == [InvalidEnumValueError(CT.single_phase_transformer, AT.tap_side, [SPT_ID], BranchSide)]


def test_single_phase_transformer_tap_pos_batch():
    spt = _spt(**{AT.tap_min: -2, AT.tap_max: 2})
    update = initialize_array(DatasetType.update, CT.single_phase_transformer, (3, 1))
    update[AT.id] = SPT_ID
    update[AT.tap_pos] = [[1], [9], [-2]]
    errors = validate_batch_data(_spt_input(spt), {CT.single_phase_transformer: update}, symmetric=False)
    assert errors is not None
    assert set(errors) == {1}
    assert errors[1] == _tap_error()


@pytest.mark.parametrize(
    ("overrides", "u_rated_to", "match"),
    [
        pytest.param({AT.from_connection: WindingConnection1P.a}, None, "u1 of single_phase_transformer", id="u1-a"),
        pytest.param({AT.u2: U_MV}, None, "u2 of single_phase_transformer", id="u2-a"),
        pytest.param(
            {AT.u2: 240.0, AT.to_connection: WindingConnection1P.ab_center_tapped},
            U_SPLIT_PHASE,
            "u2 of single_phase_transformer",
            id="u2-center-tapped",
        ),
        pytest.param({AT.u2: U_SWER}, U_MV / 2, "u2 of single_phase_transformer", id="u2-node"),
    ],
)
def test_single_phase_transformer_winding_voltage_warning(overrides: dict, u_rated_to: float | None, match: str):
    spt = _spt(**overrides)
    input_data = _spt_input(spt, u_rated_to=u_rated_to)
    with pytest.warns(ValidationWarning, match=match) as record:
        assert validate_input_data(input_data, symmetric=False) is None
    assert len(record) == 1
    assert str(SPT_ID) in str(record[0].message)
    with pytest.warns(ValidationWarning, match=match):
        assert_valid_input_data(input_data, symmetric=False)


@pytest.mark.parametrize(
    ("overrides", "u_rated_to"),
    [
        pytest.param({}, None, id="ab-and-a"),
        pytest.param({AT.u2: U_MV, AT.to_connection: WindingConnection1P.ac}, None, id="ac"),
        pytest.param({AT.u2: 480.0, AT.to_connection: WindingConnection1P.ab_center_tapped}, U_SPLIT_PHASE, id="ct"),
        pytest.param({AT.u2: 1.1 * U_SWER}, None, id="ten-percent-off"),
        pytest.param({AT.u2: 0.85 * U_SWER}, None, id="fifteen-percent-off"),
        pytest.param({AT.to_node: 99}, None, id="unknown-node"),
        pytest.param({AT.to_connection: 9}, None, id="unknown-connection"),
    ],
)
def test_single_phase_transformer_winding_voltage_no_warning(overrides: dict, u_rated_to: float | None):
    spt = _spt(**overrides)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ValidationWarning)
        errors = validate_input_data(_spt_input(spt, u_rated_to=u_rated_to), symmetric=False)
    if AT.to_node in overrides:
        assert errors == [InvalidIdError(CT.single_phase_transformer, AT.to_node, [SPT_ID], [CT.node])]
    elif AT.to_connection in overrides and overrides[AT.to_connection] == 9:  # noqa: PLR2004
        assert errors == [
            InvalidEnumValueError(CT.single_phase_transformer, AT.to_connection, [SPT_ID], WindingConnection1P)
        ]
    else:
        assert errors is None


def test_single_phase_transformer_winding_voltage_warning_once_per_batch():
    spt = _spt(**{AT.u2: U_MV})
    input_data = _spt_input(spt)
    update = initialize_array(DatasetType.update, CT.single_phase_transformer, (3, 1))
    update[AT.id] = SPT_ID
    update[AT.to_status] = [[1], [0], [1]]
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert validate_batch_data(input_data, {CT.single_phase_transformer: update}, symmetric=False) is None
    assert [warning.category for warning in record] == [ValidationWarning]
    with pytest.warns(ValidationWarning):
        assert_valid_batch_data(input_data, {CT.single_phase_transformer: update}, symmetric=False)


def _sensor(sensor_type: CT, terminal_type: MeasuredTerminalType) -> np.ndarray:
    sensor = initialize_array(DatasetType.input, sensor_type, 1)
    sensor[AT.id] = 30
    sensor[AT.measured_object] = SPT_ID
    sensor[AT.measured_terminal_type] = terminal_type
    if sensor_type in (CT.sym_power_sensor, CT.asym_power_sensor):
        sensor[AT.power_sigma] = 1e3
        sensor[AT.p_measured] = 1e3
        sensor[AT.q_measured] = 1e2
    else:
        sensor[AT.angle_measurement_type] = AngleMeasurementType.local_angle
        sensor[AT.i_sigma] = 1.0
        sensor[AT.i_angle_sigma] = 0.1
        sensor[AT.i_measured] = 10.0
        sensor[AT.i_angle_measured] = 0.0
    return sensor


SENSOR_TYPES = [CT.sym_power_sensor, CT.asym_power_sensor, CT.sym_current_sensor, CT.asym_current_sensor]


@pytest.mark.parametrize("terminal_type", [MeasuredTerminalType.branch_from, MeasuredTerminalType.branch_to])
@pytest.mark.parametrize("sensor_type", SENSOR_TYPES)
def test_single_phase_transformer_sensor(sensor_type: CT, terminal_type: MeasuredTerminalType):
    input_data = _spt_input(_spt(), components={sensor_type: _sensor(sensor_type, terminal_type)})
    assert validate_input_data(input_data, CalculationType.state_estimation, symmetric=False) is None


@pytest.mark.parametrize("sensor_type", [CT.sym_power_sensor, CT.asym_power_sensor])
def test_single_phase_transformer_sensor_wrong_terminal_type(sensor_type: CT):
    input_data = _spt_input(_spt(), components={sensor_type: _sensor(sensor_type, MeasuredTerminalType.load)})
    errors = validate_input_data(input_data, CalculationType.state_estimation, symmetric=False)
    assert errors is not None
    assert (
        InvalidIdError(
            sensor_type,
            AT.measured_object,
            [30],
            [CT.sym_load, CT.asym_load],
            {AT.measured_terminal_type: MeasuredTerminalType.load},
        )
        in errors
    )


# Tap regulator on a single_phase_transformer
REGULATOR_SPT_ID = 40
SPT_TAPS = {AT.tap_side: BranchSide.from_side, AT.tap_min: -4, AT.tap_max: 4, AT.tap_size: 275.0}


def _regulated_spt_input(
    spt: np.ndarray | None = None, u_rated_to: float | None = None, regulator: dict | None = None
) -> dict[CT, np.ndarray]:
    """
    A source at the from node, a load at the to node and a transformer_tap_regulator (to side, u_set 12.7 kV) of the
    single_phase_transformer, which has from-side taps -4..4 of 275 V.
    """
    spt = _spt(**SPT_TAPS) if spt is None else spt
    data = _spt_input(spt, u_rated_to=u_rated_to)

    source = initialize_array(DatasetType.input, CT.source, 1)
    source_values: dict[AT, float] = {AT.id: 2, AT.node: FROM_NODE, AT.status: 1, AT.u_ref: 1.0, AT.sk: 1e8}
    for source_field, source_value in source_values.items():
        source[source_field] = source_value
    load = initialize_array(DatasetType.input, CT.asym_load, 1)
    load_values: dict[AT, Any] = {
        AT.id: 3,
        AT.node: TO_NODE,
        AT.status: 1,
        AT.type: LoadGenType.const_power,
        AT.p_specified: [[10e3, 0.0, 0.0]],
        AT.q_specified: [[1e3, 0.0, 0.0]],
    }
    for load_field, load_value in load_values.items():
        load[load_field] = load_value

    regulator_array = initialize_array(DatasetType.input, CT.transformer_tap_regulator, 1)
    regulator_values: dict[AT, Any] = {
        AT.id: REGULATOR_SPT_ID,
        AT.regulated_object: SPT_ID,
        AT.status: 1,
        AT.control_side: BranchSide.to_side,
        AT.u_set: U_SWER,
        AT.u_band: 400.0,
        **(regulator or {}),
    }
    for regulator_field, regulator_value in regulator_values.items():
        regulator_array[regulator_field] = regulator_value
    return {**data, CT.source: source, CT.asym_load: load, CT.transformer_tap_regulator: regulator_array}


def _validate_regulated(data: dict[CT, np.ndarray]):
    with warnings.catch_warnings():
        warnings.simplefilter("error", ValidationWarning)
        return validate_input_data(data, CalculationType.power_flow, symmetric=False)


def test_single_phase_transformer_tap_regulator():
    u_set, u_band = 13.1e3, 300.0
    data = _regulated_spt_input(regulator={AT.u_set: u_set, AT.u_band: u_band})
    assert _validate_regulated(data) is None
    assert validate_input_data(data, symmetric=False) is None

    # mirror: the C++ core regulates it
    model = PowerGridModel(data)
    result = model._calculate_power_flow(  # noqa: SLF001
        symmetric=False,
        tap_changing_strategy=TapChangingStrategy.any_valid_tap,
        experimental_features="enabled",
    )
    u_a = result[CT.node][AT.u][1][0]
    assert abs(u_a - u_set) <= u_band / 2
    assert result[CT.transformer_tap_regulator][AT.tap_pos][0] != 0


@pytest.mark.parametrize("control_side", [BranchSide.from_side, BranchSide.to_side])
def test_single_phase_transformer_tap_regulator_control_side_valid(control_side: BranchSide):
    u_set = U_SWER if control_side == BranchSide.to_side else U_MV
    assert _validate_regulated(_regulated_spt_input(regulator={AT.control_side: control_side, AT.u_set: u_set})) is None


def test_single_phase_transformer_tap_regulator_invalid_control_side():
    # a three-winding transformer side
    data = _regulated_spt_input(regulator={AT.control_side: Branch3Side.side_3})
    assert _validate_regulated(data) == [
        InvalidAssociatedEnumValueError(
            CT.transformer_tap_regulator, [AT.control_side, AT.regulated_object], [REGULATOR_SPT_ID], [BranchSide]
        )
    ]
    # not a side at all
    errors = _validate_regulated(_regulated_spt_input(regulator={AT.control_side: 5}))
    assert errors is not None
    assert (
        InvalidEnumValueError(
            CT.transformer_tap_regulator, AT.control_side, [REGULATOR_SPT_ID], [BranchSide, Branch3Side]
        )
        in errors
    )
    # mirror
    with pytest.raises(MissingCaseForEnumError):
        PowerGridModel(data)._calculate_power_flow(  # noqa: SLF001
            symmetric=False, tap_changing_strategy="any_valid_tap", experimental_features="enabled"
        )


@pytest.mark.parametrize("regulator_status", [1, 0])
@pytest.mark.parametrize("field", [AT.tap_min, AT.tap_max, AT.tap_size])
def test_single_phase_transformer_tap_regulator_missing_tap_field(field: AT, regulator_status: int):
    # whatever the status of the regulator
    taps = {key: value for key, value in SPT_TAPS.items() if key != field}
    data = _regulated_spt_input(_spt(**taps), regulator={AT.status: regulator_status})
    assert _validate_regulated(data) == [MissingValueError(CT.single_phase_transformer, field, [SPT_ID])]

    # mirror: tap_min and tap_max are rejected by the C++ core; it cannot tell an omitted tap_size from 0
    if field == AT.tap_size:
        PowerGridModel(data)
    else:
        with pytest.raises(AutomaticTapInputError, match=rf"single_phase_transformer {SPT_ID} .* tap_min and tap_max"):
            PowerGridModel(data)


def test_single_phase_transformer_tap_regulator_missing_all_tap_fields():
    data = _regulated_spt_input(_spt())
    assert _validate_regulated(data) == [
        MissingValueError(CT.single_phase_transformer, field, [SPT_ID])
        for field in (AT.tap_min, AT.tap_max, AT.tap_size)
    ]
    with pytest.raises(AutomaticTapInputError):
        PowerGridModel(data)


def test_single_phase_transformer_tap_fields_optional_without_regulator():
    data = _regulated_spt_input(_spt())
    del data[CT.transformer_tap_regulator]
    assert _validate_regulated(data) is None
    PowerGridModel(data)

    # another single_phase_transformer is regulated: only that one needs the tap fields
    other = _spt(**{AT.id: SPT_ID + 1, AT.from_node: FROM_NODE, AT.to_node: TO_NODE})
    data = _regulated_spt_input()
    data[CT.single_phase_transformer] = np.concatenate([data[CT.single_phase_transformer], other])
    assert _validate_regulated(data) is None


def _u_set_warned(data: dict[CT, np.ndarray]) -> list:
    """The u_set warnings (a phase-to-phase to side of this network has no earth reference)"""
    return [warning for warning in _warned(data, CalculationType.power_flow) if "u_set" in str(warning.message)]


@pytest.mark.parametrize(
    ("spt_overrides", "regulator_overrides", "u_rated_to"),
    [
        pytest.param({}, {AT.u_set: U_MV}, None, id="phase-to-earth-line-voltage"),
        pytest.param({}, {AT.u_set: 1.25 * U_SWER}, None, id="phase-to-earth-high"),
        pytest.param({}, {AT.u_set: 0.75 * U_SWER}, None, id="phase-to-earth-low"),
        pytest.param({}, {AT.u_set: U_MV, AT.status: 0}, None, id="status-off"),
        pytest.param(
            {AT.to_connection: WindingConnection1P.b}, {AT.u_set: U_MV}, None, id="phase-to-earth-other-phase"
        ),
        pytest.param(
            {AT.u2: U_MV, AT.to_connection: WindingConnection1P.bc},
            {AT.u_set: U_SWER},
            None,
            id="phase-to-phase-to-side",
        ),
        pytest.param(
            {}, {AT.control_side: BranchSide.from_side, AT.u_set: U_SWER}, None, id="phase-to-phase-from-side"
        ),
        pytest.param(
            {AT.u2: 480.0, AT.to_connection: WindingConnection1P.ab_center_tapped},
            {AT.u_set: 480.0},
            U_SPLIT_PHASE,
            id="center-tapped-full-winding",
        ),
        pytest.param(
            {AT.u2: 480.0, AT.to_connection: WindingConnection1P.ab_center_tapped},
            {AT.u_set: 0.0},
            U_SPLIT_PHASE,
            id="center-tapped-zero",
        ),
    ],
)
def test_single_phase_transformer_tap_regulator_u_set_warning(
    spt_overrides: dict, regulator_overrides: dict, u_rated_to: float | None
):
    data = _regulated_spt_input(_spt(**SPT_TAPS, **spt_overrides), u_rated_to, regulator_overrides)
    record = _u_set_warned(data)
    assert len(record) == 1
    message = str(record[0].message)
    assert f"[{REGULATOR_SPT_ID}]" in message
    assert "u_set" in message
    # the data is still valid (the phase-to-phase to side of this network has no earth reference either)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ValidationWarning)
        assert validate_input_data(data, CalculationType.power_flow, symmetric=False) is None
        assert_valid_input_data(data, CalculationType.power_flow, symmetric=False)


@pytest.mark.parametrize(
    ("spt_overrides", "regulator_overrides", "u_rated_to"),
    [
        pytest.param({}, {}, None, id="phase-to-earth"),
        pytest.param({}, {AT.u_set: 1.15 * U_SWER}, None, id="phase-to-earth-fifteen-percent-high"),
        pytest.param({}, {AT.u_set: 0.85 * U_SWER}, None, id="phase-to-earth-fifteen-percent-low"),
        pytest.param({}, {AT.control_side: BranchSide.from_side, AT.u_set: U_MV}, None, id="phase-to-phase-from-side"),
        pytest.param(
            {AT.u2: U_MV, AT.to_connection: WindingConnection1P.bc}, {AT.u_set: U_MV}, None, id="phase-to-phase-to-side"
        ),
        pytest.param(
            {AT.u2: 480.0, AT.to_connection: WindingConnection1P.ab_center_tapped},
            {AT.u_set: 240.0},
            U_SPLIT_PHASE,
            id="center-tapped-leg",
        ),
        pytest.param(
            {AT.u2: 480.0, AT.to_connection: WindingConnection1P.ab_center_tapped},
            {AT.u_set: 250.0},
            U_SPLIT_PHASE,
            id="center-tapped-leg-250",
        ),
        pytest.param({}, {AT.regulated_object: 99, AT.u_set: U_MV}, None, id="unknown-regulated-object"),
        pytest.param({}, {AT.control_side: 5, AT.u_set: U_MV}, None, id="unknown-control-side"),
        pytest.param({AT.to_node: 99}, {AT.u_set: U_MV}, None, id="unknown-node"),
        pytest.param({}, {AT.u_set: np.nan}, None, id="missing-u-set"),
    ],
)
def test_single_phase_transformer_tap_regulator_u_set_no_warning(
    spt_overrides: dict, regulator_overrides: dict, u_rated_to: float | None
):
    data = _regulated_spt_input(_spt(**SPT_TAPS, **spt_overrides), u_rated_to, regulator_overrides)
    assert _u_set_warned(data) == []


def test_single_phase_transformer_tap_regulator_u_set_warning_assertion():
    data = _regulated_spt_input(regulator={AT.u_set: U_MV})
    with pytest.warns(ValidationWarning, match="u_set") as record:
        assert_valid_input_data(data, CalculationType.power_flow, symmetric=False)
    assert len(record) == 1


def test_single_phase_transformer_tap_regulator_u_set_only_checked_for_power_flow():
    data = _regulated_spt_input(regulator={AT.u_set: U_MV})
    assert len(_warned(data)) == 1
    assert len(_warned(data, CalculationType.power_flow)) == 1
    assert _warned(data, CalculationType.state_estimation) == []
    assert _warned(data, CalculationType.short_circuit) == []


def _u_set_update(u_set: list[float]) -> dict[CT, np.ndarray]:
    update = initialize_array(DatasetType.update, CT.transformer_tap_regulator, (len(u_set), 1))
    update[AT.id] = REGULATOR_SPT_ID
    update[AT.u_set] = [[value] for value in u_set]
    return {CT.transformer_tap_regulator: update}


def _u_set_batch_warnings(data: dict[CT, np.ndarray], update: dict[CT, np.ndarray]) -> list:
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert validate_batch_data(data, update, CalculationType.power_flow, symmetric=False) is None
    return [warning for warning in record if issubclass(warning.category, ValidationWarning)]


def test_single_phase_transformer_tap_regulator_u_set_warning_batch():
    # the input data is nominal, and only scenario 1 updates u_set far off
    data = _regulated_spt_input()
    record = _u_set_batch_warnings(data, _u_set_update([U_SWER, U_MV, 1.1 * U_SWER]))
    assert len(record) == 1
    message = str(record[0].message)
    assert f"[{REGULATOR_SPT_ID}]" in message
    assert "1 of 3" in message
    assert "scenario 1" in message
    with pytest.warns(ValidationWarning, match="u_set"):
        assert_valid_batch_data(
            data, _u_set_update([U_SWER, U_MV, 1.1 * U_SWER]), CalculationType.power_flow, symmetric=False
        )

    # one warning for several scenarios
    record = _u_set_batch_warnings(data, _u_set_update([U_MV, U_SWER, 2 * U_SWER]))
    assert len(record) == 1
    assert "2 of 3" in str(record[0].message)
    assert "scenario 0" in str(record[0].message)

    # all scenarios nominal
    assert _u_set_batch_warnings(data, _u_set_update([U_SWER, 1.1 * U_SWER, 0.9 * U_SWER])) == []

    # the input data is off, but every scenario updates u_set
    off_data = _regulated_spt_input(regulator={AT.u_set: U_MV})
    assert _u_set_batch_warnings(off_data, _u_set_update([U_SWER, U_SWER])) == []


# Load and generator connection
LOAD_GEN_ID = 5
REGULATOR_ID = 6
ASYM_LOAD_GEN_TYPES = [CT.asym_load, CT.asym_gen]
SYM_LOAD_GEN_TYPES = [CT.sym_load, CT.sym_gen]
LOAD_GEN_TYPES = ASYM_LOAD_GEN_TYPES + SYM_LOAD_GEN_TYPES
NA_INT_S = np.iinfo(np.int8).min
INVALID_CONNECTIONS = [2, 5]


def _load_gen(component: CT, connection: int, load_gen_type: LoadGenType = LoadGenType.const_power) -> np.ndarray:
    load_gen = initialize_array(DatasetType.input, component, 1)
    load_gen[AT.id] = LOAD_GEN_ID
    load_gen[AT.node] = 2
    load_gen[AT.status] = 1
    load_gen[AT.type] = load_gen_type
    load_gen[AT.connection] = connection
    if component in (CT.sym_load, CT.sym_gen):
        load_gen[AT.p_specified] = 1e3
        load_gen[AT.q_specified] = 1e2
    else:
        load_gen[AT.p_specified] = [[1e3, 1e3, 1e3]]
        load_gen[AT.q_specified] = [[1e2, 1e2, 1e2]]
    return load_gen


def _load_gen_input(load_gen: np.ndarray, component: CT, regulated: bool = False) -> dict[CT, np.ndarray]:
    node = initialize_array(DatasetType.input, CT.node, 2)
    node[AT.id] = [1, 2]
    node[AT.u_rated] = 10e3
    line = initialize_array(DatasetType.input, CT.line, 1)
    line[AT.id] = 3
    line[AT.from_node] = 1
    line[AT.to_node] = 2
    line[AT.from_status] = 1
    line[AT.to_status] = 1
    line[AT.r1] = 0.1
    line[AT.x1] = 0.1
    line[AT.c1] = 1e-8
    line[AT.tan1] = 0.0
    line[AT.r0] = 0.3
    line[AT.x0] = 0.3
    line[AT.c0] = 1e-8
    line[AT.tan0] = 0.0
    source = initialize_array(DatasetType.input, CT.source, 1)
    source[AT.id] = 4
    source[AT.node] = 1
    source[AT.status] = 1
    source[AT.u_ref] = 1.0
    data = {CT.node: node, CT.line: line, CT.source: source, component: load_gen}
    if regulated:
        regulator = initialize_array(DatasetType.input, CT.voltage_regulator, 1)
        regulator[AT.id] = REGULATOR_ID
        regulator[AT.regulated_object] = LOAD_GEN_ID
        regulator[AT.status] = 1
        regulator[AT.u_ref] = 1.0
        data[CT.voltage_regulator] = regulator
    return data


def _validate_load_gen(component: CT, connection: int, **kwargs):
    load_gen = _load_gen(component, connection, **kwargs)
    return validate_input_data(_load_gen_input(load_gen, component), symmetric=False)


def _regulation_error() -> InvalidVoltageRegulationError:
    return InvalidVoltageRegulationError(CT.voltage_regulator, AT.regulated_object, [REGULATOR_ID])


@pytest.mark.parametrize("connection", [LoadGenConnection.wye, LoadGenConnection.delta, NA_INT_S])
@pytest.mark.parametrize("component", ASYM_LOAD_GEN_TYPES)
def test_asym_load_gen_valid_connection(component: CT, connection: int):
    assert _validate_load_gen(component, connection) is None


@pytest.mark.parametrize("connection", INVALID_CONNECTIONS)
@pytest.mark.parametrize("component", LOAD_GEN_TYPES)
def test_load_gen_invalid_connection(component: CT, connection: int):
    assert _validate_load_gen(component, connection) == [
        InvalidEnumValueError(component, AT.connection, [LOAD_GEN_ID], LoadGenConnection)
    ]


@pytest.mark.parametrize("component", SYM_LOAD_GEN_TYPES)
def test_sym_load_gen_wye_connection(component: CT):
    assert _validate_load_gen(component, LoadGenConnection.wye) is None


@pytest.mark.parametrize("component", SYM_LOAD_GEN_TYPES)
def test_sym_load_gen_delta_connection(component: CT):
    load_gen = _load_gen(component, LoadGenConnection.delta)
    input_data = _load_gen_input(load_gen, component)
    assert validate_input_data(input_data, symmetric=False) == [
        InvalidValueError(component, AT.connection, [LOAD_GEN_ID], [LoadGenConnection.wye])
    ]
    with pytest.raises(InvalidPhaseConnection):
        PowerGridModel(input_data)


@pytest.mark.parametrize("load_gen_type", [LoadGenType.const_power, LoadGenType.const_impedance])
@pytest.mark.parametrize("component", ASYM_LOAD_GEN_TYPES)
def test_voltage_regulator_on_delta_asym_load_gen(component: CT, load_gen_type: LoadGenType):
    load_gen = _load_gen(component, LoadGenConnection.delta, load_gen_type)
    input_data = _load_gen_input(load_gen, component, regulated=True)
    errors = validate_input_data(input_data, CalculationType.power_flow, symmetric=False)
    # exactly one error for the regulator, also if the load_gen is not const_power either
    assert errors == [_regulation_error()]


@pytest.mark.parametrize("component", ASYM_LOAD_GEN_TYPES)
def test_voltage_regulator_on_delta_asym_load_gen_mirror(component: CT):
    load_gen = _load_gen(component, LoadGenConnection.delta)
    input_data = _load_gen_input(load_gen, component, regulated=True)
    model = PowerGridModel(input_data)
    with pytest.raises(UnsupportedLoadGenTypeForVoltageRegulator):
        model._calculate_power_flow(symmetric=False, experimental_features="enabled")  # noqa: SLF001


@pytest.mark.parametrize("component", ASYM_LOAD_GEN_TYPES)
def test_voltage_regulator_on_wye_asym_load_gen(component: CT):
    load_gen = _load_gen(component, LoadGenConnection.wye)
    input_data = _load_gen_input(load_gen, component, regulated=True)
    assert validate_input_data(input_data, CalculationType.power_flow, symmetric=False) is None
    model = PowerGridModel(input_data)
    model._calculate_power_flow(symmetric=False, experimental_features="enabled")  # noqa: SLF001


# phase-to-phase windings without an earth reference (D67)
FLOAT_SPT_ID = 151
FLOAT_LOAD_ID = 301
FLOAT_LINE_ID = 401
FLOAT_TRANSFORMER_ID = 501
FLOAT_NO_REFERENCE = "no earth reference"
NODE_PRIMARY, NODE_SOURCE, NODE_SECONDARY, NODE_EXTRA = 0, 1, 3, 4


def _floating_spt(**overrides) -> np.ndarray:
    """An 11 kV / 415 V phase-to-phase transformer from the source node (1) to the secondary node (3)"""
    return _spt(
        **{
            AT.id: FLOAT_SPT_ID,
            AT.from_node: NODE_SOURCE,
            AT.to_node: NODE_SECONDARY,
            AT.u1: 11e3,
            AT.u2: 415.0,
            AT.sn: 50e3,
            AT.uk: 0.04,
            AT.pk: 600.0,
            AT.to_connection: WindingConnection1P.ab,
            **overrides,
        }
    )


def _floating_load(load_type: LoadGenType, component: CT = CT.asym_load, **overrides) -> np.ndarray:
    load = initialize_array(DatasetType.input, component, 1)
    values = {
        AT.id: FLOAT_LOAD_ID,
        AT.node: NODE_SECONDARY,
        AT.status: 1,
        AT.type: load_type,
        AT.p_specified: [[20e3, 12e3, 0.0]],
        AT.q_specified: [[5e3, 3e3, 0.0]],
        **overrides,
    }
    for field, value in values.items():
        load[field] = value
    return load


def _three_phase_transformer(winding_from: WindingType, winding_to: WindingType, clock: int) -> np.ndarray:
    """A 11 kV / 11 kV transformer from the source node (1) to the primary node (0)"""
    transformer = initialize_array(DatasetType.input, CT.transformer, 1)
    values = {
        AT.id: FLOAT_TRANSFORMER_ID,
        AT.from_node: NODE_SOURCE,
        AT.to_node: NODE_PRIMARY,
        AT.from_status: 1,
        AT.to_status: 1,
        AT.u1: 11e3,
        AT.u2: 11e3,
        AT.sn: 1e6,
        AT.uk: 0.05,
        AT.pk: 5e3,
        AT.i0: 0.001,
        AT.p0: 500.0,
        AT.winding_from: winding_from,
        AT.winding_to: winding_to,
        AT.clock: clock,
        AT.tap_side: BranchSide.from_side,
        AT.tap_pos: 0,
        AT.tap_min: -1,
        AT.tap_max: 1,
        AT.tap_size: 0.0,
    }
    for field, value in values.items():
        transformer[field] = value
    return transformer


def _floating_shunt(node: int) -> np.ndarray:
    shunt = initialize_array(DatasetType.input, CT.shunt, 1)
    for field, value in {
        AT.id: 601,
        AT.node: node,
        AT.status: 1,
        AT.g1: 1e-6,
        AT.b1: 1e-6,
        AT.g0: 1e-6,
        AT.b0: 1e-6,
    }.items():
        shunt[field] = value
    return shunt


def _floating_input(
    load_type: LoadGenType = LoadGenType.const_power,
    *,
    load_component: CT = CT.asym_load,
    load_overrides: dict | None = None,
    spt_overrides: dict | None = None,
    extra: dict[CT, np.ndarray] | None = None,
    with_load: bool = True,
    source_status: int | None = 1,
) -> dict[CT, np.ndarray]:
    """
    A source at 11 kV (node 1), a phase-to-phase transformer (ab to ab) to the 415 V node 3 and an asym_load there.
    Node 0 (11 kV) and node 4 (415 V) are free for the variants. source_status None means no source.
    """
    node = initialize_array(DatasetType.input, CT.node, 4)
    node[AT.id] = [NODE_PRIMARY, NODE_SOURCE, NODE_SECONDARY, NODE_EXTRA]
    node[AT.u_rated] = [11e3, 11e3, 415.0, 415.0]
    data: dict[CT, np.ndarray] = {CT.node: node, CT.single_phase_transformer: _floating_spt(**(spt_overrides or {}))}
    if source_status is not None:
        source = initialize_array(DatasetType.input, CT.source, 1)
        for field, value in {
            AT.id: 2,
            AT.node: NODE_SOURCE,
            AT.status: source_status,
            AT.u_ref: 1.0,
            AT.sk: 1e8,
        }.items():
            source[field] = value
        data[CT.source] = source
    if with_load:
        data[load_component] = _floating_load(load_type, load_component, **(load_overrides or {}))
    for component, array in (extra or {}).items():
        data[component] = np.concatenate([data[component], array]) if component in data else array
    return data


def _warned(data: dict[CT, np.ndarray], calculation_type: CalculationType | None = None) -> list:
    """The ValidationWarnings that validate_input_data emits for this data"""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        validate_input_data(data, calculation_type, symmetric=False)
    return [warning for warning in record if issubclass(warning.category, ValidationWarning)]


def test_floating_winding_warns():
    data = _floating_input()
    with pytest.warns(ValidationWarning, match=FLOAT_NO_REFERENCE) as record:
        assert validate_input_data(data, symmetric=False) is None
    assert len(record) == 1
    message = str(record[0].message)
    assert str(FLOAT_SPT_ID) in message
    assert "to side" in message
    assert "from side" not in message
    with pytest.warns(ValidationWarning, match=FLOAT_NO_REFERENCE):
        assert_valid_input_data(data, symmetric=False)


def _asym_line_from_secondary() -> np.ndarray:
    asym_line = _swer_asym_line(PhaseConnection.ab)
    asym_line[AT.id] = FLOAT_LINE_ID
    asym_line[AT.from_node] = NODE_SECONDARY
    asym_line[AT.to_node] = NODE_EXTRA
    return asym_line


@pytest.mark.parametrize(
    ("kwargs", "calculation_type", "side"),
    [
        pytest.param({"load_type": LoadGenType.const_current}, None, "to", id="const-current"),
        pytest.param(
            {"load_type": LoadGenType.const_impedance, "load_overrides": {AT.p_specified: 0.0, AT.q_specified: 0.0}},
            None,
            "to",
            id="const-impedance-zero-power",
        ),
        pytest.param(
            {"load_type": LoadGenType.const_impedance, "load_overrides": {AT.status: 0}}, None, "to", id="status-off"
        ),
        pytest.param(
            {"load_type": LoadGenType.const_impedance, "load_overrides": {AT.connection: LoadGenConnection.delta}},
            None,
            "to",
            id="delta",
        ),
        pytest.param(
            {"load_type": LoadGenType.const_power, "extra": {CT.asym_line: _asym_line_from_secondary()}},
            None,
            "to",
            id="line-capacitance",
        ),
        pytest.param(
            {"load_type": LoadGenType.const_impedance}, CalculationType.short_circuit, "to", id="short-circuit"
        ),
        pytest.param(
            {
                "load_type": LoadGenType.const_impedance,
                "spt_overrides": {AT.from_node: NODE_PRIMARY, AT.to_connection: WindingConnection1P.a, AT.u2: 240.0},
                "extra": {CT.transformer: _three_phase_transformer(WindingType.wye_n, WindingType.delta, 1)},
            },
            None,
            "from",
            id="primary-behind-delta",
        ),
    ],
)
def test_floating_winding_not_anchored_by(kwargs: dict, calculation_type: CalculationType | None, side: str):
    record = _warned(_floating_input(**kwargs), calculation_type)
    assert len(record) == 1
    assert FLOAT_NO_REFERENCE in str(record[0].message)
    assert f"{FLOAT_SPT_ID} ({side} side)" in str(record[0].message)


def _primary_variant(winding_from: WindingType, winding_to: WindingType, clock: int, **kwargs) -> dict:
    """Node 0 is fed through a three-phase transformer and transformer 151 goes from node 0 (ab) to node 3 (a)"""
    return {
        "load_type": LoadGenType.const_power,
        "spt_overrides": {AT.from_node: NODE_PRIMARY, AT.to_connection: WindingConnection1P.a, AT.u2: 240.0},
        "extra": {CT.transformer: _three_phase_transformer(winding_from, winding_to, clock), **kwargs.pop("extra", {})},
        **kwargs,
    }


@pytest.mark.parametrize(
    "kwargs",
    [
        pytest.param(
            {"load_type": LoadGenType.const_impedance, "load_overrides": {AT.connection: LoadGenConnection.wye}},
            id="const-impedance-asym-load",
        ),
        pytest.param(
            {"load_type": LoadGenType.const_impedance, "load_component": CT.asym_gen}, id="const-impedance-asym-gen"
        ),
        pytest.param(
            {
                "extra": {
                    CT.single_phase_transformer: _floating_spt(
                        **{
                            AT.id: FLOAT_SPT_ID + 1,
                            AT.from_node: NODE_SECONDARY,
                            AT.to_node: NODE_EXTRA,
                            AT.from_connection: WindingConnection1P.a,
                            AT.to_connection: WindingConnection1P.a,
                            AT.u1: 415.0 / np.sqrt(3.0),
                            AT.u2: 415.0 / np.sqrt(3.0),
                        }
                    )
                }
            },
            id="phase-to-earth-winding",
        ),
        pytest.param(
            {
                "spt_overrides": {AT.to_connection: WindingConnection1P.ab_center_tapped, AT.u2: 830.0 / np.sqrt(3.0)},
                "with_load": False,
            },
            id="center-tapped",
        ),
        pytest.param(_primary_variant(WindingType.wye_n, WindingType.wye_n, 0), id="primary-wye-n"),
        pytest.param(_primary_variant(WindingType.delta, WindingType.wye_n, 11), id="primary-delta-wye-n"),
        pytest.param(_primary_variant(WindingType.delta, WindingType.zigzag_n, 0), id="primary-zigzag-n"),
        pytest.param(
            _primary_variant(WindingType.wye_n, WindingType.delta, 1, extra={CT.shunt: _floating_shunt(NODE_PRIMARY)}),
            id="primary-shunt",
        ),
    ],
)
def test_floating_winding_anchored_by(kwargs: dict):
    data = _floating_input(**kwargs)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ValidationWarning)
        assert validate_input_data(data, symmetric=False) is None


@pytest.mark.parametrize(
    "kwargs",
    [
        pytest.param({"source_status": None}, id="no-source"),
        pytest.param({"source_status": 0}, id="source-off"),
        pytest.param({"spt_overrides": {AT.from_status: 0}}, id="from-side-off"),
    ],
)
def test_floating_winding_not_energised(kwargs: dict):
    assert _warned(_floating_input(**kwargs)) == []


def _load_status_update(statuses: list[list[int]]) -> dict[CT, np.ndarray]:
    update = initialize_array(DatasetType.update, CT.asym_load, (len(statuses), 1))
    update[AT.id] = FLOAT_LOAD_ID
    update[AT.status] = statuses
    return {CT.asym_load: update}


def _batch_warnings(data: dict[CT, np.ndarray], update: dict[CT, np.ndarray]) -> list:
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert validate_batch_data(data, update, symmetric=False) is None
    return [warning for warning in record if issubclass(warning.category, ValidationWarning)]


def test_floating_winding_batch():
    # the input data is anchored, and only scenario 1 switches the load off
    data = _floating_input(LoadGenType.const_impedance)
    update = _load_status_update([[1], [0], [1]])
    record = _batch_warnings(data, update)
    assert len(record) == 1
    message = str(record[0].message)
    assert FLOAT_NO_REFERENCE in message
    assert "1 of 3" in message
    assert "scenario 1" in message
    with pytest.warns(ValidationWarning, match=FLOAT_NO_REFERENCE):
        assert_valid_batch_data(data, update, symmetric=False)

    # all scenarios anchored
    assert _batch_warnings(data, _load_status_update([[1], [1], [1]])) == []

    # the input data is floating, but every scenario switches the load on
    floating_input = _floating_input(LoadGenType.const_impedance, load_overrides={AT.status: 0})
    assert _batch_warnings(floating_input, _load_status_update([[1], [1], [1]])) == []

    # one warning for several scenarios
    record = _batch_warnings(data, _load_status_update([[0], [1], [0]]))
    assert len(record) == 1
    assert "2 of 3" in str(record[0].message)
    assert "scenario 0" in str(record[0].message)


def _calculate_floating_power_flow(model: PowerGridModel, calculation_method: CalculationMethod):
    return model._calculate_power_flow(  # noqa: SLF001
        symmetric=False, calculation_method=calculation_method, experimental_features="enabled"
    )


def test_floating_winding_mirror():
    # the warned network: the iterative current method does not converge
    data = _floating_input()
    assert len(_warned(data)) == 1
    with pytest.raises(IterationDiverge):
        _calculate_floating_power_flow(PowerGridModel(data), CalculationMethod.iterative_current)

    # with a const_impedance load there is no warning and the methods agree
    data = _floating_input(LoadGenType.const_impedance)
    assert _warned(data) == []
    model = PowerGridModel(data)
    iterative = _calculate_floating_power_flow(model, CalculationMethod.iterative_current)
    newton = _calculate_floating_power_flow(model, CalculationMethod.newton_raphson)
    np.testing.assert_allclose(iterative[CT.node][AT.u], newton[CT.node][AT.u], rtol=1e-6)


def _reference_components(n_node: int, from_idx: np.ndarray, to_idx: np.ndarray) -> np.ndarray:
    """A plain BFS reference: the smallest node index of the component of every node"""
    neighbours: list[list[int]] = [[] for _ in range(n_node)]
    for first, second in zip(from_idx.tolist(), to_idx.tolist(), strict=True):
        neighbours[first].append(second)
        neighbours[second].append(first)
    labels = [-1] * n_node
    for start in range(n_node):
        if labels[start] >= 0:
            continue
        labels[start] = start
        queue = [start]
        while queue:
            current = queue.pop()
            for other in neighbours[current]:
                if labels[other] < 0:
                    labels[other] = start
                    queue.append(other)
    return np.array(labels)


def _edges(*pairs: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    return np.array([pair[0] for pair in pairs], dtype=np.intp), np.array([pair[1] for pair in pairs], dtype=np.intp)


@pytest.mark.parametrize(
    ("n_node", "edges", "expected"),
    [
        pytest.param(5, [(0, 1), (1, 2), (2, 3), (3, 4)], [0, 0, 0, 0, 0], id="chain"),
        pytest.param(5, [(4, 3), (3, 2), (2, 1), (1, 0)], [0, 0, 0, 0, 0], id="reversed-chain"),
        pytest.param(6, [(0, 1), (3, 4), (4, 5)], [0, 0, 2, 3, 3, 3], id="two-islands"),
        pytest.param(4, [], [0, 1, 2, 3], id="no-edges"),
        pytest.param(4, [(1, 2)], [0, 1, 1, 3], id="isolated-nodes"),
        pytest.param(3, [(1, 1), (2, 2)], [0, 1, 2], id="self-loop"),
        pytest.param(0, [], [], id="no-nodes"),
    ],
)
def test_connected_components(n_node: int, edges: list[tuple[int, int]], expected: list[int]):
    from_idx, to_idx = _edges(*edges)
    np.testing.assert_array_equal(_connected_components(n_node, from_idx, to_idx), expected)


@pytest.mark.parametrize("seed", range(5))
def test_connected_components_random(seed: int):
    rng = np.random.default_rng(seed)
    n_node, n_edge = 50, 40
    from_idx = rng.integers(0, n_node, n_edge)
    to_idx = rng.integers(0, n_node, n_edge)
    labels = _connected_components(n_node, from_idx, to_idx)
    reference = _reference_components(n_node, from_idx, to_idx)
    # the same partition, and every label is the smallest index of its component
    _, inverse = np.unique(labels, return_inverse=True)
    _, reference_inverse = np.unique(reference, return_inverse=True)
    np.testing.assert_array_equal(inverse, reference_inverse)
    np.testing.assert_array_equal(labels, reference)
