# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""
The Python validation accepts every reduced-phase, split-phase and delta validation case in tests/data.

The cases are selected by their content: a single_phase_transformer, an asym_line with phases, or a load or generator
with a connection. New cases are picked up automatically.

The expected-error cases (e.g. delta-on-absent-phase, swer-fault-on-absent-phase and missing-phase-for-winding) are
included on purpose: their errors depend on the topology (the phases present at each node), which only the C++ core
checks. The Python validation checks each component's own data, so it must accept them. The exception is
missing-tap-limits: a regulated single_phase_transformer without both tap limits is caught by the Python validation as
well (D73), so that case must give exactly that error.

No case may raise a ValidationWarning either: the winding voltages of every single_phase_transformer are nominal, and
every phase-to-phase winding has an earth reference (D67).
"""

import json
import warnings
from pathlib import Path

import pytest

from power_grid_model import AttributeType, CalculationType, ComponentType
from power_grid_model.utils import json_deserialize_from_file
from power_grid_model.validation import ValidationWarning, validate_batch_data, validate_input_data
from power_grid_model.validation.errors import MissingValueError

DATA_PATH = Path(__file__).parents[2] / "data"
MIN_CASE_COUNT = 53
# the cases whose errors the Python validation reports as well, with the expected errors
EXPECTED_ERRORS = {
    "power_flow/single-phase-transformer-tap-regulator/missing-tap-limits": [
        MissingValueError(ComponentType.single_phase_transformer, AttributeType.tap_max, [154])
    ],
}


def _reduced_phase_cases():
    for path in sorted(DATA_PATH.glob("*/**/input.json")):
        raw = json.loads(path.read_text(encoding="utf-8"))["data"]
        if (
            "single_phase_transformer" in raw
            or any("phases" in row for row in raw.get("asym_line", []))
            or any(
                "connection" in row
                for component in ("sym_load", "sym_gen", "asym_load", "asym_gen")
                for row in raw.get(component, [])
            )
        ):
            yield pytest.param(path.parent, id=str(path.parent.relative_to(DATA_PATH)).replace("\\", "/"))


CASES = list(_reduced_phase_cases())


def test_case_count():
    assert len(CASES) >= MIN_CASE_COUNT
    case_ids = {case.id for case in CASES}
    assert set(EXPECTED_ERRORS) <= case_ids


@pytest.mark.parametrize("case_dir", CASES)
def test_reduced_phase_case_is_valid(case_dir: Path):
    calculation_type = CalculationType[case_dir.relative_to(DATA_PATH).parts[0]]
    expected_errors = EXPECTED_ERRORS.get(case_dir.relative_to(DATA_PATH).as_posix())
    input_data = json_deserialize_from_file(case_dir / "input.json")
    with warnings.catch_warnings():
        warnings.simplefilter("error", ValidationWarning)
        assert validate_input_data(input_data, calculation_type, symmetric=False) == expected_errors
        if (case_dir / "update_batch.json").exists():
            update_data = json_deserialize_from_file(case_dir / "update_batch.json")
            assert validate_batch_data(input_data, update_data, calculation_type, symmetric=False) is None
