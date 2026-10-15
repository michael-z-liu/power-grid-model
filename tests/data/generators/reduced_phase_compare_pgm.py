# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""Compare PGM with the OpenDSS references of the reduced-phase validation cases.

For every case and method it prints, per field, the max absolute and relative deviation, and the worst ratio
|actual - expected| / (atol + rtol |expected|) for the tolerances in params.json and for candidate tolerances.
The fields are transformed exactly like tests/unit/utils.py compare_result does (angles as phasors, u_angle
normalised on the first entry). A ratio below 0.1 means the tolerance passes with a 10x margin.

The root must contain only cases that this report understands, so generate the cases into a scratch directory first
and compare there (tests/data itself also holds other validation cases). Run from the repository root with the
installed package:

    uv run --no-project --with opendssdirect.py --with numpy \\
        python tests/data/generators/reduced_phase_opendss_cases.py <scratch dir>
    uv run --no-sync python tests/data/generators/reduced_phase_compare_pgm.py <scratch dir>
"""

import json
import re
import sys
from pathlib import Path

import numpy as np

from power_grid_model import PowerGridModel
from power_grid_model.enum import CalculationMethod, ShortCircuitVoltageScaling, _ExperimentalFeatures
from power_grid_model.utils import json_deserialize_from_file

CANDIDATES = [(r, a) for r in (1e-10, 1e-9, 1e-8, 1e-7, 1e-6) for a in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4)]


def fields(actual, reference):
    """(component.field, actual, expected) like compare_result"""
    for component, rows in reference.items():
        act = actual[component]
        for name in rows[0]:
            if name == "id":
                continue
            exp = np.array([row[name] for row in rows], dtype=float)
            a = np.asarray(act[name], dtype=float)
            if name == "u_angle":
                a = a - a.ravel()[0]
                exp = exp - exp.ravel()[0]
            if name.endswith("_angle"):
                mag = name[: -len("_angle")]
                a = np.asarray(act[mag]) * np.exp(1j * a)
                exp = np.array([row[mag] for row in rows], dtype=float) * np.exp(1j * exp)
            yield f"{component}.{name}", a, exp


def atol_for(atol, name):
    if not isinstance(atol, dict):
        return atol
    for pattern, value in atol.items():
        if re.fullmatch(pattern, name.split(".")[1]):
            return value
    return atol["default"]


def worst_ratio(items, rtol, atol):
    return max(float(np.max(np.abs(a - e) / (atol_for(atol, name) + rtol * np.abs(e)))) for name, a, e in items)


def run(case, params, method):
    model = PowerGridModel(json_deserialize_from_file(case / "input.json"))
    enabled = _ExperimentalFeatures.enabled
    if (case / "asym_output.json").exists():
        result = model._calculate_power_flow(
            symmetric=False, calculation_method=CalculationMethod[method], experimental_features=enabled
        )
        return [(result, json.loads((case / "asym_output.json").read_text())["data"])]
    scaling = ShortCircuitVoltageScaling[params["short_circuit_voltage_scaling"]]
    if (case / "sc_output.json").exists():
        result = model._calculate_short_circuit(short_circuit_voltage_scaling=scaling, experimental_features=enabled)
        return [(result, json.loads((case / "sc_output.json").read_text())["data"])]
    update = json_deserialize_from_file(case / "update_batch.json")
    result = model._calculate_short_circuit(
        update_data=update, short_circuit_voltage_scaling=scaling, experimental_features=enabled
    )
    references = json.loads((case / "sc_output_batch.json").read_text())["data"]
    return [({k: v[s] for k, v in result.items()}, ref) for s, ref in enumerate(references)]


def main():
    root = Path(sys.argv[1])
    for params_path in sorted(root.rglob("params.json")):
        case = params_path.parent
        params = json.loads(params_path.read_text())
        if "raises" in params:
            continue
        methods = params["calculation_method"]
        methods = methods if isinstance(methods, list) else [methods]
        print(case.relative_to(root), f"(params: rtol {params['rtol']}, atol {params['atol']})")
        all_items = []
        for method in methods:
            for k, (actual, reference) in enumerate(run(case, params, method)):
                items = list(fields(actual, reference))
                all_items += items
                parts = []
                for name, a, e in items:
                    dev = np.abs(a - e)
                    rel = np.where(np.abs(e) > 1e-3, dev / np.maximum(np.abs(e), 1e-300), 0.0)
                    parts.append(f"{name} {dev.max():.1e}/{rel.max():.1e}")
                ratio = worst_ratio(items, params["rtol"], params["atol"])
                print(f"  {method} #{k}: ratio {ratio:.3f} | abs/rel: " + ", ".join(parts))
        print(f"  overall ratio with params: {worst_ratio(all_items, params['rtol'], params['atol']):.3f}")
        passing = [(r, a) for r, a in CANDIDATES if worst_ratio(all_items, r, a) < 0.1]
        print("  candidates with 10x margin (rtol, atol):", ", ".join(f"({r:.0e}, {a:.0e})" for r, a in passing[:8]))


if __name__ == "__main__":
    main()
