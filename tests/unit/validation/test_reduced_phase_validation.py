# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""
Validation of reduced-phase (one- or two-phase) networks.

The "mirror" tests check that an input flagged by the Python validation is also rejected by the C++ core, so that the
two do not drift apart.
"""

import warnings

import numpy as np
import pytest

from power_grid_model import (
    AttributeType as AT,
    ComponentType as CT,
    DatasetType,
    PowerGridModel,
    initialize_array,
)
from power_grid_model.enum import (
    AngleMeasurementType,
    BranchSide,
    CalculationType,
    MeasuredTerminalType,
    PhaseConnection,
    WindingConnection1P,
)
from power_grid_model.errors import (
    InvalidPhaseConnection,
    InvalidTransformerClock,
    MissingCaseForEnumError,
)
from power_grid_model.validation import (
    ValidationWarning,
    assert_valid_batch_data,
    assert_valid_input_data,
    validate_batch_data,
    validate_input_data,
)
from power_grid_model.validation.errors import (
    InvalidEnumValueError,
    InvalidIdError,
    InvalidValueError,
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
    errors = _validate_spt(_spt(**{field: value}))
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


def test_single_phase_transformer_tap_regulator():
    regulator = initialize_array(DatasetType.input, CT.transformer_tap_regulator, 1)
    regulator[AT.id] = 40
    regulator[AT.regulated_object] = SPT_ID
    regulator[AT.status] = 1
    regulator[AT.control_side] = BranchSide.to_side
    regulator[AT.u_set] = U_SWER
    regulator[AT.u_band] = 200.0
    input_data = _spt_input(_spt(), components={CT.transformer_tap_regulator: regulator})
    errors = validate_input_data(input_data, CalculationType.power_flow, symmetric=False)
    assert errors is not None
    assert (
        InvalidIdError(
            CT.transformer_tap_regulator,
            AT.regulated_object,
            [40],
            [CT.transformer, CT.three_winding_transformer],
        )
        in errors
    )
