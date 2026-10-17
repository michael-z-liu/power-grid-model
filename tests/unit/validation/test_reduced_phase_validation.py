# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""
Validation of reduced-phase (one- or two-phase) networks.

The "mirror" tests check that an input flagged by the Python validation is also rejected by the C++ core, so that the
two do not drift apart.
"""

import numpy as np
import pytest

from power_grid_model import (
    AttributeType as AT,
    ComponentType as CT,
    DatasetType,
    PowerGridModel,
    initialize_array,
)
from power_grid_model.enum import PhaseConnection
from power_grid_model.errors import InvalidPhaseConnection, MissingCaseForEnumError
from power_grid_model.validation import validate_batch_data, validate_input_data
from power_grid_model.validation.errors import (
    InvalidEnumValueError,
    MissingValueError,
    MultiFieldValidationError,
    NotApplicableValueError,
    NotGreaterThanError,
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
