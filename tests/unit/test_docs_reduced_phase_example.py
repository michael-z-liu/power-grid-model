# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""
Run the worked example of docs/user_manual/reduced-phase-networks.md.

The example is the Python code between the markers <!-- example-start --> and <!-- example-end --> of that page. It
uses the private _calculate_power_flow, so this test also fails when that signature changes.
"""

import re
import warnings
from pathlib import Path

import numpy as np
import pytest

from power_grid_model import ComponentType
from power_grid_model.errors import ExperimentalFeature
from power_grid_model.validation import ValidationWarning

DOC = Path(__file__).parents[2] / "docs" / "user_manual" / "reduced-phase-networks.md"
START_MARKER = "<!-- example-start -->"
END_MARKER = "<!-- example-end -->"

LEG_VOLTAGE_RANGE = (230.0, 245.0)
SWER_ANGLE = 30.0
ANGLE_TOLERANCE = 1.0


def _example_code() -> str:
    text = DOC.read_text(encoding="utf-8")
    assert text.count(START_MARKER) == 1
    assert text.count(END_MARKER) == 1
    section = text.split(START_MARKER)[1].split(END_MARKER)[0]
    blocks = re.findall(r"^```python\n(.*?)^```$", section, flags=re.DOTALL | re.MULTILINE)
    assert len(blocks) > 1
    return "\n".join(blocks)


@pytest.fixture(scope="module")
def example() -> dict:
    namespace: dict = {}
    with warnings.catch_warnings():
        warnings.simplefilter("error", ValidationWarning)
        exec(compile(_example_code(), str(DOC), "exec"), namespace)
    return namespace


def test_example_input_is_valid(example: dict):
    assert example["errors"] is None


def test_example_is_gated_without_experimental_features(example: dict):
    assert isinstance(example["gated"], ExperimentalFeature)


def test_example_absent_phases_are_zero(example: dict):
    u = example["u"]
    assert u.shape == (4, 3)
    np.testing.assert_array_equal(u[1:3, 1:], 0.0)  # SWER nodes 2 and 3: only phase a
    assert u[3, 2] == 0.0  # split-phase node 4: no phase c


def test_example_split_phase_legs(example: dict):
    u = example["u"]
    u_angle = example["u_angle"]
    low, high = LEG_VOLTAGE_RANGE
    assert low < u[3, 0] < high
    assert low < u[3, 1] < high
    # leg b is opposite to leg a
    assert abs((u_angle[3, 0] - u_angle[3, 1]) % 360.0 - 180.0) < ANGLE_TOLERANCE


def test_example_swer_angle(example: dict):
    u_angle = example["u_angle"]
    assert u_angle[1, 0] == pytest.approx(SWER_ANGLE, abs=ANGLE_TOLERANCE)
    assert u_angle[2, 0] == pytest.approx(SWER_ANGLE, abs=ANGLE_TOLERANCE)


def test_example_delta_load_is_split_over_the_legs(example: dict):
    load = example["result"][ComponentType.asym_load]
    delta_p = load["p"][1]
    assert delta_p[0] == pytest.approx(2e3, rel=1e-2)
    assert delta_p[1] == pytest.approx(2e3, rel=1e-2)
    assert delta_p[2] == 0.0
