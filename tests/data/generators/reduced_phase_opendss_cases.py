# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""Generator for the reduced-phase asym_line validation cases, with OpenDSS as the reference solver.

Each case is defined once. The script writes the PGM input (input.json, params.json, update_batch.json), the OpenDSS
model (test.dss) and the reference output computed by OpenDSS, plus a .license file per file. It generates:

- tests/data/power_flow/asym-line/{single-phase-a, two-phase-ab, swer-a, reduced-phase-feeder,
  reduced-phase-feeder-linear, reduced-phase-symmetric-not-supported}
- tests/data/state_estimation/reduced-phase-not-supported
- tests/data/short_circuit/{asym-line-reduced-phase, asym-line-swer-a}

Run from the repository root (opendssdirect.py is not a project dependency, so it is added on the fly):

    uv run --no-project --with opendssdirect.py --with numpy \\
        python tests/data/generators/reduced_phase_opendss_cases.py tests/data

Optionally pass case path fragments after the root to generate only the matching cases.

Modelling choices:

- source: Vsource with Z1/Z0 in ohm from sk, rx_ratio and z01_ratio, exactly as PGM computes them; in a short
  circuit the voltage is c * u_rated (pu=1.1) and the impedance is c * u_rated^2 / sk (IEC 60909, as in PGM)
- lines: explicit R/X/C matrices of the present phases; a line with a neutral gets the numpy Kron-reduced matrix,
  which is cross-checked against an OpenDSS LineCode with Kron=yes
- loads: every PGM load (also a sym_load) becomes single-phase wye loads with kV = u_rated / sqrt(3) exactly;
  const_pq -> model=1, const_y -> model=2, const_i -> model=5; a generator is a load with negative power
- a branch side with status 0 is connected to its own dummy bus in OpenDSS (the open-end pi model); PGM reports
  the isolated node and the open side as 0, so the reference writes 0 for them as well
"""

import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import opendssdirect as dss

FREQ = 50.0
# REUSE-IgnoreStart
LICENSE = (
    "SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>\n"
    "\n"
    "SPDX-License-Identifier: MPL-2.0\n"
)
# REUSE-IgnoreEnd
PHASE_CONNECTION = {"abc": 0, "a": 1, "b": 2, "c": 3, "ab": 4, "ac": 5, "bc": 6}
LOAD_TYPE = {"const_pq": 0, "const_y": 1, "const_i": 2}
DSS_LOAD_MODEL = {"const_pq": 1, "const_y": 2, "const_i": 5}
FAULT_TYPE = {"three_phase": 0, "single_phase_to_ground": 1, "two_phase": 2, "two_phase_to_ground": 3}
FAULT_PHASE = {"abc": 0, "a": 1, "b": 2, "c": 3, "ab": 4, "ac": 5, "bc": 6}
CONDUCTORS = "abcn"

# ---------------------------------------------------------------------------------------------------------------
# line data per km


def matrix4(self_aa, self_nn, mutual_pp, mutual_pn):
    """symmetric 4x4 (a, b, c, n) with the given self and mutual values; mutual_pp = (ab, ac, bc)"""
    m = np.zeros((4, 4))
    m[0, 0] = m[1, 1] = m[2, 2] = self_aa
    m[3, 3] = self_nn
    m[0, 1] = m[1, 0] = mutual_pp[0]
    m[0, 2] = m[2, 0] = mutual_pp[1]
    m[1, 2] = m[2, 1] = mutual_pp[2]
    for p in range(3):
        m[p, 3] = m[3, p] = mutual_pn[p]
    return m


R_PER_KM = matrix4(0.25, 0.35, (0.050, 0.048, 0.050), (0.049, 0.047, 0.051))
X_PER_KM = matrix4(0.40, 0.45, (0.150, 0.130, 0.150), (0.140, 0.120, 0.145))
C_PER_KM = np.array([[9.0e-9, -1.5e-9, -1.2e-9], [-1.5e-9, 9.2e-9, -1.5e-9], [-1.2e-9, -1.5e-9, 9.0e-9]])
C0_PER_KM = 5.0e-9
C1_PER_KM = 10.0e-9

# SWER line, see the design doc: Carson earth return with rho = 100 ohm m, GMR = 1.27 mm
R_SWER_PER_KM = 1.40 + math.pi**2 * FREQ * 1e-4
X_SWER_PER_KM = 4.0 * math.pi * FREQ * 1e-4 * math.log(658.5 * math.sqrt(100.0 / FREQ) / 0.00127)
C_SWER_PER_KM = 6.0e-9


@dataclass
class Line:
    id: int
    from_node: int
    to_node: int
    phases: str  # present phases, e.g. "a", "bc", "abc"
    length: float  # km
    neutral: bool = False
    c_form: str = "matrix"  # "matrix" or "c0c1" (abc only)
    swer: bool = False
    from_status: int = 1
    to_status: int = 1

    def idx(self):
        return ["abc".index(p) for p in self.phases]

    def conductor_idx(self):
        return self.idx() + ([3] if self.neutral else [])

    def z_full(self):
        """series impedance (ohm) of the conductors (present phases and the neutral)"""
        if self.swer:
            return np.array([[complex(R_SWER_PER_KM, X_SWER_PER_KM)]]) * self.length
        ci = self.conductor_idx()
        return (R_PER_KM + 1j * X_PER_KM)[np.ix_(ci, ci)] * self.length

    def z_phase(self):
        """Kron-reduced series impedance (ohm) of the present phases"""
        z = self.z_full()
        if not self.neutral:
            return z
        n = len(self.phases)
        return z[:n, :n] - np.outer(z[:n, n], z[n, :n]) / z[n, n]

    def c_phase(self):
        """capacitance (F) of the present phases"""
        if self.swer:
            return np.array([[C_SWER_PER_KM]]) * self.length
        if self.c_form == "c0c1":
            assert self.phases == "abc"
            cs = (2.0 * C1_PER_KM + C0_PER_KM) / 3.0
            cm = (C0_PER_KM - C1_PER_KM) / 3.0
            return (np.full((3, 3), cm) + np.eye(3) * (cs - cm)) * self.length
        i = self.idx()
        return C_PER_KM[np.ix_(i, i)] * self.length

    def pgm(self):
        result = {
            "id": self.id,
            "from_node": self.from_node,
            "to_node": self.to_node,
            "from_status": self.from_status,
            "to_status": self.to_status,
        }
        if self.phases != "abc":
            result["phases"] = PHASE_CONNECTION[self.phases]
        z = self.z_full()
        ci = self.conductor_idx()
        for row in range(4):
            for col in range(row + 1):
                name = CONDUCTORS[row] + CONDUCTORS[col]
                if row in ci and col in ci:
                    value = z[ci.index(row), ci.index(col)]
                    result["r_" + name] = value.real
                    result["x_" + name] = value.imag
        if self.c_form == "c0c1":
            result["c0"] = C0_PER_KM * self.length
            result["c1"] = C1_PER_KM * self.length
        else:
            c = self.c_phase()
            pi_ = self.idx()
            for row in range(3):
                for col in range(row + 1):
                    if row in pi_ and col in pi_:
                        result["c_" + CONDUCTORS[row] + CONDUCTORS[col]] = c[pi_.index(row), pi_.index(col)]
        return result


@dataclass
class Load:
    id: int
    node: int
    kind: str  # asym_load, asym_gen, sym_load
    type: str
    p: list  # per phase (W) for asym, total for sym
    q: list
    status: int = 1

    def pgm(self):
        result = {"id": self.id, "node": self.node, "status": self.status, "type": LOAD_TYPE[self.type]}
        if self.kind == "sym_load":
            result["p_specified"] = self.p
            result["q_specified"] = self.q
        else:
            result["p_specified"] = list(self.p)
            result["q_specified"] = list(self.q)
        return result

    def per_phase(self):
        """injection-free per-phase consumption (W, var) as seen by OpenDSS"""
        if self.kind == "sym_load":
            return [self.p / 3.0] * 3, [self.q / 3.0] * 3
        sign = -1.0 if self.kind == "asym_gen" else 1.0
        return [sign * x for x in self.p], [sign * x for x in self.q]


@dataclass
class Fault:
    id: int
    node: int
    fault_type: str
    fault_phase: str
    r_f: float
    x_f: float = 0.0


@dataclass
class Case:
    path: str  # relative to tests/data
    u_rated: float
    nodes: list
    lines: list
    loads: list = field(default_factory=list)
    calculation: str = "power_flow"  # power_flow or short_circuit
    params: dict = field(default_factory=dict)
    faults: list = field(default_factory=list)  # short circuit: one fault per scenario, batch if more than one
    # cases that only check an error: (output type, raises name)
    raises_only: tuple | None = None
    source_id: int = 201
    sk: float = 100e6
    rx_ratio: float = 0.1
    z01_ratio: float = 1.0


# ---------------------------------------------------------------------------------------------------------------
# PGM input


def pgm_input(case: Case, fault: Fault | None):
    data = {
        "node": [{"id": n, "u_rated": case.u_rated} for n in case.nodes],
        "asym_line": [line.pgm() for line in case.lines],
        "source": [
            {
                "id": case.source_id,
                "node": case.nodes[0],
                "status": 1,
                "u_ref": 1.0,
                "sk": case.sk,
                "rx_ratio": case.rx_ratio,
                "z01_ratio": case.z01_ratio,
            }
        ],
    }
    for kind in ("sym_load", "asym_load", "asym_gen"):
        items = [load.pgm() for load in case.loads if load.kind == kind]
        if items:
            data[kind] = items
    if fault is not None:
        data["fault"] = [pgm_fault(fault, with_status=True)]
    return {"version": "1.0", "type": "input", "is_batch": False, "attributes": {}, "data": data}


def pgm_fault(fault: Fault, with_status: bool):
    result: dict[str, int | float] = {"id": fault.id}
    if with_status:
        result["status"] = 1
    result.update(
        {
            "fault_type": FAULT_TYPE[fault.fault_type],
            "fault_phase": FAULT_PHASE[fault.fault_phase],
            "fault_object": fault.node,
            "r_f": fault.r_f,
            "x_f": fault.x_f,
        }
    )
    return result


# ---------------------------------------------------------------------------------------------------------------
# OpenDSS model


def fmt_matrix(m):
    """lower triangle in OpenDSS notation"""
    rows = [" ".join(repr(float(m[row, col])) for col in range(row + 1)) for row in range(m.shape[0])]
    return "[" + " | ".join(rows) + "]"


def nodes_of(phases):
    return "".join("." + str("abc".index(p) + 1) for p in phases)


def source_z(case: Case):
    z_abs = case.u_rated**2 / case.sk
    x1 = z_abs / math.sqrt(case.rx_ratio**2 + 1.0)
    r1 = x1 * case.rx_ratio
    return complex(r1, x1), complex(r1, x1) * case.z01_ratio


def dss_script(case: Case, fault: Fault | None, voltage_pu: float, with_loads: bool):
    z1, z0 = source_z(case)
    if fault is not None:
        # IEC 60909: sk is defined at the rated voltage, so the source impedance is c * u_rated^2 / sk (as in PGM)
        z1 *= voltage_pu
        z0 *= voltage_pu
    kv = case.u_rated / 1e3
    lines = [
        "Clear",
        f"Set DefaultBaseFrequency={FREQ:g}",
        (
            f"New Circuit.c phases=3 basekv={kv!r} pu={voltage_pu!r} angle=0 bus1=n{case.nodes[0]}"
            f" Z1=[{z1.real!r}, {z1.imag!r}] Z0=[{z0.real!r}, {z0.imag!r}]"
        ),
    ]
    for line in case.lines:
        z = line.z_phase()
        c_nf = line.c_phase() * 1e9
        bus1 = f"n{line.from_node}" if line.from_status else f"open_l{line.id}_from"
        bus2 = f"n{line.to_node}" if line.to_status else f"open_l{line.id}_to"
        nodes = nodes_of(line.phases)
        lines.append(
            f"New Line.l{line.id} phases={len(line.phases)} bus1={bus1}{nodes} bus2={bus2}{nodes} length=1 units=km"
            f" rmatrix={fmt_matrix(z.real)} xmatrix={fmt_matrix(z.imag)} cmatrix={fmt_matrix(c_nf)}"
        )
    if with_loads:
        kv_phase = case.u_rated / math.sqrt(3.0) / 1e3
        for load in case.loads:
            if not load.status:
                continue
            p, q = load.per_phase()
            for k in range(3):
                if p[k] == 0.0 and q[k] == 0.0:
                    continue
                lines.append(
                    f"New Load.ld{load.id}_{'abc'[k]} bus1=n{load.node}.{k + 1} phases=1 conn=wye"
                    f" model={DSS_LOAD_MODEL[load.type]} kV={kv_phase!r} kW={p[k] / 1e3!r} kvar={q[k] / 1e3!r}"
                    " vminpu=0.1 vmaxpu=2"
                )
    if fault is not None:
        if fault.fault_type == "single_phase_to_ground":
            k = "abc".index(fault.fault_phase) + 1
            lines.append(f"New Fault.f{fault.id} phases=1 bus1=n{fault.node}.{k} r={fault.r_f!r}")
        elif fault.fault_type == "two_phase":
            k1, k2 = ("abc".index(p) + 1 for p in fault.fault_phase)
            lines.append(
                f"New Fault.f{fault.id} phases=1 bus1=n{fault.node}.{k1} bus2=n{fault.node}.{k2} r={fault.r_f!r}"
            )
        elif fault.fault_type == "three_phase":
            # three phase-to-earth resistances r_f, as PGM models a three-phase fault with impedance
            lines.append(
                f"New Fault.f{fault.id} phases=3 bus1=n{fault.node}.1.2.3 bus2=n{fault.node}.0.0.0 r={fault.r_f!r}"
            )
        else:
            raise NotImplementedError(fault.fault_type)
    lines += [f"Set VoltageBases=[{kv!r}]", "CalcVoltageBases", "Set Tolerance=1e-12", "Set MaxIterations=200"]
    if fault is not None:
        lines.append("Set Mode=Snapshot")
    lines.append("Solve")
    return lines


def run_dss(script):
    dss.Basic.ClearAll()
    for command in script:
        dss.Text.Command(command)
    assert dss.Solution.Converged(), "OpenDSS did not converge"


def check_kron_against_opendss(case: Case):
    """the numpy Kron reduction equals an OpenDSS LineCode with Kron=yes"""
    for line in case.lines:
        if not line.neutral:
            continue
        z4 = line.z_full()
        n = z4.shape[0]
        dss.Basic.ClearAll()
        dss.Text.Command("Clear")
        dss.Text.Command(f"Set DefaultBaseFrequency={FREQ:g}")
        dss.Text.Command("New Circuit.k basekv=10 bus1=x")
        dss.Text.Command(
            f"New LineCode.k{line.id} nphases={n} units=km rmatrix={fmt_matrix(z4.real)} xmatrix={fmt_matrix(z4.imag)}"
            f" cmatrix={fmt_matrix(np.zeros((n, n)))} Neutral={n} Kron=yes"
        )
        dss.LineCodes.Name(f"k{line.id}")
        assert dss.LineCodes.Phases() == n - 1, f"OpenDSS did not Kron-reduce line {line.id}"
        r = np.array(dss.LineCodes.Rmatrix()).reshape(n - 1, n - 1)
        x = np.array(dss.LineCodes.Xmatrix()).reshape(n - 1, n - 1)
        expected = line.z_phase()
        dev = np.max(np.abs(r + 1j * x - expected))
        assert dev < 1e-12 * np.max(np.abs(expected)), f"Kron mismatch for line {line.id}: {dev}"
        print(f"  Kron check line {line.id}: OpenDSS LineCode agrees with numpy (max dev {dev:.1e} ohm)")


# ---------------------------------------------------------------------------------------------------------------
# reference output


def bus_voltages(node):
    dss.Circuit.SetActiveBus(f"n{node}")
    u = [0.0, 0.0, 0.0]
    angle = [0.0, 0.0, 0.0]
    mag_angle = dss.Bus.VMagAngle()
    for k, number in enumerate(dss.Bus.Nodes()):
        if 1 <= number <= 3:
            u[number - 1] = mag_angle[2 * k]
            angle[number - 1] = math.radians(mag_angle[2 * k + 1])
    return u, angle


def terminal_values(name, values_getter, n_cond):
    """per terminal, per phase complex values of the active element"""
    dss.Circuit.SetActiveElement(name)
    raw = values_getter()
    order = dss.CktElement.NodeOrder()
    result = [np.zeros(3, dtype=complex), np.zeros(3, dtype=complex)]
    for terminal in range(2):
        for k in range(n_cond):
            idx = terminal * n_cond + k
            number = order[idx]
            if 1 <= number <= 3:
                result[terminal][number - 1] = complex(raw[2 * idx], raw[2 * idx + 1])
    return result


def isolated_nodes(case: Case):
    """nodes that are only reached through an open branch side (PGM reports them as 0)"""
    energized = {case.nodes[0]}
    changed = True
    while changed:
        changed = False
        for line in case.lines:
            if line.from_status and line.to_status:
                for a, b in ((line.from_node, line.to_node), (line.to_node, line.from_node)):
                    if a in energized and b not in energized:
                        energized.add(b)
                        changed = True
    return [n for n in case.nodes if n not in energized]


def mag_angle(values):
    return [float(abs(v)) for v in values], [float(np.angle(v)) if abs(v) > 0.0 else 0.0 for v in values]


def pf_reference(case: Case):
    isolated = isolated_nodes(case)
    nodes = []
    for node in case.nodes:
        if node in isolated:
            nodes.append({"id": node, "u": [0.0, 0.0, 0.0], "u_angle": [0.0, 0.0, 0.0]})
            continue
        u, angle = bus_voltages(node)
        nodes.append({"id": node, "u": u, "u_angle": angle})
    lines = []
    for line in case.lines:
        s_from, s_to = terminal_values(f"Line.l{line.id}", dss.CktElement.Powers, len(line.phases))
        s_from = s_from * 1e3 if line.from_status else np.zeros(3, dtype=complex)
        s_to = s_to * 1e3 if line.to_status else np.zeros(3, dtype=complex)
        lines.append(
            {
                "id": line.id,
                "p_from": [float(v.real) for v in s_from],
                "q_from": [float(v.imag) for v in s_from],
                "p_to": [float(v.real) for v in s_to],
                "q_to": [float(v.imag) for v in s_to],
            }
        )
    return {"node": nodes, "asym_line": lines}


def sc_reference(case: Case, fault: Fault):
    isolated = isolated_nodes(case)
    nodes = []
    for node in case.nodes:
        if node in isolated:
            nodes.append({"id": node, "u": [0.0, 0.0, 0.0], "u_angle": [0.0, 0.0, 0.0]})
            continue
        u, angle = bus_voltages(node)
        nodes.append({"id": node, "u": u, "u_angle": angle})
    lines = []
    for line in case.lines:
        i_from, i_to = terminal_values(f"Line.l{line.id}", dss.CktElement.Currents, len(line.phases))
        if not line.from_status:
            i_from = np.zeros(3, dtype=complex)
        if not line.to_status:
            i_to = np.zeros(3, dtype=complex)
        i_from_mag, i_from_angle = mag_angle(i_from)
        i_to_mag, i_to_angle = mag_angle(i_to)
        lines.append(
            {
                "id": line.id,
                "i_from": i_from_mag,
                "i_from_angle": i_from_angle,
                "i_to": i_to_mag,
                "i_to_angle": i_to_angle,
            }
        )
    # fault current: the current from the node into the fault on each faulted phase
    dss.Circuit.SetActiveElement(f"Fault.f{fault.id}")
    raw = dss.CktElement.Currents()
    order = dss.CktElement.NodeOrder()
    n_cond = dss.CktElement.NumConductors()
    i_f = np.zeros(3, dtype=complex)
    for terminal in range(2):
        for k in range(n_cond):
            idx = terminal * n_cond + k
            number = order[idx]
            if 1 <= number <= 3:
                value = complex(raw[2 * idx], raw[2 * idx + 1])
                if fault.fault_type == "two_phase" or terminal == 0:
                    i_f[number - 1] = value
    i_f_mag, i_f_angle = mag_angle(i_f)
    return {"node": nodes, "asym_line": lines, "fault": [{"id": fault.id, "i_f": i_f_mag, "i_f_angle": i_f_angle}]}


# ---------------------------------------------------------------------------------------------------------------
# writing


def write_json(path: Path, content):
    path.write_text(json.dumps(content, indent=2) + "\n")
    Path(str(path) + ".license").write_text(LICENSE)


def write_text(path: Path, content: str):
    path.write_text(content)
    Path(str(path) + ".license").write_text(LICENSE)


def output_dataset(output_type, data, is_batch=False):
    return {"version": "1.0", "type": output_type, "is_batch": is_batch, "attributes": {}, "data": data}


def generate(case: Case, root: Path):
    target = root / case.path
    target.mkdir(parents=True, exist_ok=True)
    print(f"{case.path}")
    first_fault = case.faults[0] if case.faults else None
    write_json(target / "input.json", pgm_input(case, first_fault))
    write_json(target / "params.json", case.params)

    if case.raises_only is not None:
        output_type, _ = case.raises_only
        component_ids = {"node": case.nodes, "asym_line": [line.id for line in case.lines]}
        write_json(
            target / f"{output_type}.json",
            output_dataset(output_type, {k: [{"id": i} for i in v] for k, v in component_ids.items()}),
        )
        return

    check_kron_against_opendss(case)

    if case.calculation == "power_flow":
        script = dss_script(case, None, 1.0, with_loads=True)
        run_dss(script)
        write_text(target / "test.dss", "\n".join(script) + "\n")
        write_json(target / "asym_output.json", output_dataset("asym_output", pf_reference(case)))
        return

    scripts = []
    references = []
    for fault in case.faults:
        script = dss_script(case, fault, 1.1, with_loads=False)
        run_dss(script)
        scripts.append(script)
        references.append(sc_reference(case, fault))
    dss_text = []
    for k, script in enumerate(scripts):
        if len(scripts) > 1:
            dss_text.append(f"! scenario {k}")
        dss_text += script
        dss_text.append("")
    write_text(target / "test.dss", "\n".join(dss_text))
    if len(case.faults) == 1:
        write_json(target / "sc_output.json", output_dataset("sc_output", references[0]))
    else:
        write_json(
            target / "update_batch.json",
            {
                "version": "1.0",
                "type": "update",
                "is_batch": True,
                "attributes": {},
                "data": [{"fault": [pgm_fault(fault, with_status=False)]} for fault in case.faults],
            },
        )
        write_json(target / "sc_output_batch.json", output_dataset("sc_output", references, is_batch=True))


# ---------------------------------------------------------------------------------------------------------------
# cases

# Tolerances: per case, the smallest round rtol (with atol <= 1e-6), then the smallest round atol, for which the worst
# |actual - expected| / (atol + rtol |expected|) is below 0.1 (a 10x margin); see compare_reduced_phase_cases.py.
# PGM and OpenDSS agree to about 3.4e-10 relative everywhere; the cases with iterative_current are limited by its
# convergence criterion (error_tolerance 1e-8 on the voltage gives up to 2e-8 relative on the power).
SC_PARAMS = {
    "calculation_method": "iec60909",
    "rtol": 1e-8,
    "atol": 1e-8,
    "short_circuit_voltage_scaling": "maximum",
    "experimental_features": "enabled",
}


def pf_params(methods, rtol, atol):
    return {"calculation_method": methods, "rtol": rtol, "atol": atol, "experimental_features": "enabled"}


def feeder_lines():
    return [
        Line(101, 1, 2, "abc", 2.0),
        Line(102, 2, 3, "abc", 3.0, neutral=True, c_form="c0c1"),
        Line(103, 2, 4, "a", 2.0),
        Line(104, 2, 5, "bc", 2.5, neutral=True),
        Line(105, 5, 6, "c", 1.5),
        Line(106, 3, 7, "ac", 2.0),
        Line(107, 3, 8, "b", 1.0),
        Line(108, 4, 9, "a", 1.0),
        Line(109, 9, 7, "a", 1.5),
        Line(110, 8, 10, "b", 0.8, to_status=0),
    ]


def feeder_loads(all_const_y):
    def t(load_type):
        return "const_y" if all_const_y else load_type

    return [
        Load(301, 4, "asym_load", t("const_pq"), [60e3, 0.0, 0.0], [20e3, 0.0, 0.0]),
        Load(302, 5, "asym_load", t("const_i"), [0.0, 50e3, 40e3], [0.0, 15e3, 10e3]),
        Load(303, 6, "asym_load", t("const_y"), [0.0, 0.0, 30e3], [0.0, 0.0, 8e3]),
        Load(304, 7, "asym_load", t("const_pq"), [40e3, 0.0, 35e3], [10e3, 0.0, 10e3]),
        Load(305, 8, "asym_load", t("const_i"), [0.0, 25e3, 0.0], [0.0, 6e3, 0.0]),
        Load(306, 9, "asym_load", t("const_y"), [30e3, 0.0, 0.0], [9e3, 0.0, 0.0]),
        Load(401, 7, "asym_gen", t("const_pq"), [20e3, 0.0, 0.0], [5e3, 0.0, 0.0]),
        Load(501, 3, "sym_load", t("const_i"), 150e3, 50e3),
    ]


FEEDER_NODES = list(range(1, 11))
SWER_U = 22e3


def swer_lines():
    return [
        Line(101, 1, 2, "a", 20.0, swer=True),
        Line(102, 2, 3, "a", 10.0, swer=True),
        Line(103, 2, 4, "a", 5.0, swer=True, to_status=0),
    ]


def small_reduced_case(path, raises_only, params):
    return Case(
        path=path,
        u_rated=10e3,
        nodes=[1, 2],
        lines=[Line(101, 1, 2, "a", 5.0)],
        loads=[Load(301, 2, "asym_load", "const_pq", [200e3, 0.0, 0.0], [50e3, 0.0, 0.0])],
        params=params,
        raises_only=raises_only,
    )


CASES = [
    Case(
        path="power_flow/asym-line/single-phase-a",
        u_rated=10e3,
        nodes=[1, 2],
        lines=[Line(101, 1, 2, "a", 5.0)],
        loads=[Load(301, 2, "asym_load", "const_pq", [200e3, 0.0, 0.0], [50e3, 0.0, 0.0])],
        params=pf_params(["newton_raphson", "iterative_current"], 1e-6, 1e-8),
    ),
    Case(
        path="power_flow/asym-line/two-phase-ab",
        u_rated=10e3,
        nodes=[1, 2],
        lines=[Line(101, 1, 2, "ab", 3.0, neutral=True)],
        loads=[Load(301, 2, "asym_load", "const_y", [150e3, 100e3, 0.0], [40e3, 30e3, 0.0])],
        params=pf_params(["linear", "linear_current", "newton_raphson", "iterative_current"], 1e-8, 1e-8),
    ),
    Case(
        path="power_flow/asym-line/swer-a",
        u_rated=SWER_U,
        nodes=[1, 2, 3, 4],
        lines=swer_lines(),
        loads=[
            Load(301, 2, "asym_load", "const_pq", [30e3, 0.0, 0.0], [10e3, 0.0, 0.0]),
            Load(302, 3, "asym_load", "const_pq", [20e3, 0.0, 0.0], [5e3, 0.0, 0.0]),
        ],
        params=pf_params(["newton_raphson", "iterative_current"], 1e-7, 1e-7),
    ),
    Case(
        path="power_flow/asym-line/reduced-phase-feeder",
        u_rated=10e3,
        nodes=FEEDER_NODES,
        lines=feeder_lines(),
        loads=feeder_loads(all_const_y=False),
        params=pf_params(["newton_raphson", "iterative_current"], 1e-6, 1e-6),
    ),
    Case(
        path="power_flow/asym-line/reduced-phase-feeder-linear",
        u_rated=10e3,
        nodes=FEEDER_NODES,
        lines=feeder_lines(),
        loads=feeder_loads(all_const_y=True),
        params=pf_params(["linear", "linear_current"], 1e-8, 1e-8),
    ),
    small_reduced_case(
        "power_flow/asym-line/reduced-phase-symmetric-not-supported",
        ("sym_output", "ReducedPhaseNotSupported"),
        {
            "calculation_method": "newton_raphson",
            "rtol": 1e-6,
            "atol": 1e-4,
            "experimental_features": "enabled",
            "raises": {
                "raises": "ReducedPhaseNotSupported",
                "reason": "Symmetric power flow is not supported in a network with one- or two-phase nodes.",
            },
        },
    ),
    small_reduced_case(
        "state_estimation/reduced-phase-not-supported",
        ("asym_output", "ReducedPhaseNotSupported"),
        {
            "calculation_method": ["iterative_linear", "newton_raphson"],
            "rtol": 1e-6,
            "atol": 1e-4,
            "experimental_features": "enabled",
            "raises": {
                "raises": "ReducedPhaseNotSupported",
                "reason": "State estimation is not supported in a network with one- or two-phase nodes.",
            },
        },
    ),
    Case(
        path="short_circuit/asym-line-reduced-phase",
        u_rated=10e3,
        nodes=FEEDER_NODES,
        lines=feeder_lines(),
        calculation="short_circuit",
        params=SC_PARAMS,
        faults=[
            Fault(601, 9, "single_phase_to_ground", "a", 0.5),
            Fault(601, 6, "single_phase_to_ground", "c", 0.5),
            Fault(601, 5, "two_phase", "bc", 0.5),
            Fault(601, 3, "three_phase", "abc", 0.5),
        ],
    ),
    Case(
        path="short_circuit/asym-line-swer-a",
        u_rated=SWER_U,
        nodes=[1, 2, 3, 4],
        lines=swer_lines(),
        calculation="short_circuit",
        params=SC_PARAMS,
        faults=[Fault(601, 3, "single_phase_to_ground", "a", 10.0)],
    ),
]


def main():
    root = Path(sys.argv[1])
    only = sys.argv[2:]
    for case in CASES:
        if only and not any(o in case.path for o in only):
            continue
        generate(case, root)


if __name__ == "__main__":
    main()
