# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""Compare PGM with the OpenDSS references of the reduced-phase validation cases.

For every case and method it prints, per field, the max absolute and relative deviation, and the worst ratio
|actual - expected| / (atol + rtol |expected|) for the tolerances in params.json and for candidate tolerances.
The fields are transformed exactly like tests/unit/utils.py compare_result does (angles as phasors, u_angle
normalised on the first entry). A ratio below 0.1 means the tolerance passes with a 10x margin.

Per calculation method the parameters are those of params.json overridden by params["extra_params"][<method>], as in
the validation harnesses: rtol, atol, error_tolerance (default 1e-8), max_iterations (default 20) and
calculation_initialization (default "default"). The "ratio with params" of a method uses its own rtol and atol. The
candidate search is over all methods together, so it is not affected by extra_params.

If no candidate of the standard grid (rtol up to 1e-6, atol up to 1e-4) passes, it also prints the passing
candidates of a wider grid (rtol up to 1e-5, atol up to 1e-2). Reference rows may be in the compact format (lists, with
the root "attributes"). It handles single power flows and short circuits, and batches of both (update_batch.json with
asym_output_batch.json or sc_output_batch.json).

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
WIDER_CANDIDATES = [
    (r, a) for r in (1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5) for a in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
]


def reference_data(path):
    """the data of a reference dataset (a list of scenarios for a batch), with compact rows (lists) mapped to dicts
    by the root "attributes\""""
    raw = json.loads(path.read_text())
    attributes = raw.get("attributes", {})

    def scenario(data):
        return {
            component: [dict(zip(attributes[component], row)) if isinstance(row, list) else row for row in rows]
            for component, rows in data.items()
        }

    return [scenario(data) for data in raw["data"]] if raw["is_batch"] else scenario(raw["data"])


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


def method_params(params, method):
    """the params of one calculation method: the top-level params overridden by extra_params[<method>]"""
    return dict(params, **params.get("extra_params", {}).get(method, {}))


def run(case, params, method):
    mp = method_params(params, method)
    solver_options = {
        "error_tolerance": mp.get("error_tolerance", 1e-8),
        "max_iterations": mp.get("max_iterations", 20),
        "calculation_initialization": mp.get("calculation_initialization", "default"),
    }
    model = PowerGridModel(json_deserialize_from_file(case / "input.json"))
    enabled = _ExperimentalFeatures.enabled
    if (case / "asym_output.json").exists():
        result = model._calculate_power_flow(
            symmetric=False,
            calculation_method=CalculationMethod[method],
            **solver_options,
            experimental_features=enabled,
        )
        return [(result, reference_data(case / "asym_output.json"))]
    if (case / "asym_output_batch.json").exists():
        update = json_deserialize_from_file(case / "update_batch.json")
        result = model._calculate_power_flow(
            symmetric=False,
            calculation_method=CalculationMethod[method],
            **solver_options,
            update_data=update,
            experimental_features=enabled,
        )
        references = reference_data(case / "asym_output_batch.json")
        return [({k: v[s] for k, v in result.items()}, ref) for s, ref in enumerate(references)]
    scaling = ShortCircuitVoltageScaling[params["short_circuit_voltage_scaling"]]
    if (case / "sc_output.json").exists():
        result = model._calculate_short_circuit(short_circuit_voltage_scaling=scaling, experimental_features=enabled)
        return [(result, reference_data(case / "sc_output.json"))]
    update = json_deserialize_from_file(case / "update_batch.json")
    result = model._calculate_short_circuit(
        update_data=update, short_circuit_voltage_scaling=scaling, experimental_features=enabled
    )
    references = reference_data(case / "sc_output_batch.json")
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
        overall = 0.0
        for method in methods:
            for k, (actual, reference) in enumerate(run(case, params, method)):
                items = list(fields(actual, reference))
                all_items += items
                parts = []
                for name, a, e in items:
                    dev = np.abs(a - e)
                    rel = np.where(np.abs(e) > 1e-3, dev / np.maximum(np.abs(e), 1e-300), 0.0)
                    parts.append(f"{name} {dev.max():.1e}/{rel.max():.1e}")
                mp = method_params(params, method)
                ratio = worst_ratio(items, mp["rtol"], mp["atol"])
                overall = max(overall, ratio)
                print(f"  {method} #{k}: ratio {ratio:.3f} | abs/rel: " + ", ".join(parts))
        print(f"  overall ratio with params: {overall:.3f}")
        passing = [(r, a) for r, a in CANDIDATES if worst_ratio(all_items, r, a) < 0.1]
        print("  candidates with 10x margin (rtol, atol):", ", ".join(f"({r:.0e}, {a:.0e})" for r, a in passing[:8]))
        if not passing:
            wider = [(r, a) for r, a in WIDER_CANDIDATES if worst_ratio(all_items, r, a) < 0.1]
            print(
                "  wider candidates with 10x margin (rtol, atol):",
                ", ".join(f"({r:.0e}, {a:.0e})" for r, a in wider[:8]) or "none",
            )


if __name__ == "__main__":
    main()
