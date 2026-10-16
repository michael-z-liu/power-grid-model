# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""Generator for the reduced-phase asym_line validation cases, with OpenDSS as the reference solver.

Each case is defined once. The script writes the PGM input (input.json, params.json, update_batch.json), the OpenDSS
model (test.dss) and the reference output computed by OpenDSS, plus a .license file per file. It generates:

- tests/data/power_flow/asym-line/{single-phase-a, two-phase-ab, two-phase-bc-neutral-opendss-kron, swer-a,
  reduced-phase-feeder, reduced-phase-feeder-linear, reduced-phase-symmetric-not-supported}
- tests/data/state_estimation/reduced-phase-not-supported
- tests/data/short_circuit/{asym-line-reduced-phase, asym-line-swer-a}
- tests/data/power_flow/single-phase-transformer/{swer-chain, phase-to-phase, phase-to-phase-linear,
  missing-phase-for-winding}
- tests/data/short_circuit/single-phase-transformer-swer

Run from the repository root (opendssdirect.py is not a project dependency, so it is added on the fly):

    uv run --no-project --with opendssdirect.py --with numpy \\
        python tests/data/generators/reduced_phase_opendss_cases.py tests/data

Optionally pass case path fragments after the root to generate only the matching cases.

Modelling choices:

- source: Vsource with Z1/Z0 in ohm from sk, rx_ratio and z01_ratio, exactly as PGM computes them; in a short
  circuit the voltage is c * u_rated (pu=1.1) and the impedance is c * u_rated^2 / sk (IEC 60909, as in PGM)
- faults: OpenDSS replaces a fault resistance of 0 by 1e-4 ohm, so a bolted fault (r_f = 0) gets 1e-13 ohm instead,
  which is below the OpenDSS agreement of about 3e-10 relative
- lines: explicit R/X/C matrices of the present phases; a line with a neutral gets the numpy Kron-reduced matrix,
  which is cross-checked against an OpenDSS LineCode with Kron=yes; the exception is a line with opendss_kron=True,
  which gets the full matrices including the neutral as a LineCode with Kron=yes, so OpenDSS reduces the neutral itself
- loads: every PGM load (also a sym_load) becomes single-phase wye loads with kV = u_rated / sqrt(3) exactly;
  const_pq -> model=1, const_y -> model=2, const_i -> model=5; a generator is a load with negative power
- a branch side with status 0 is connected to its own dummy bus in OpenDSS (the open-end pi model); PGM reports
  the isolated node and the open side as 0, so the reference writes 0 for them as well
- single-phase transformers: a 1-phase, 2-winding wye-wye Transformer whose windings are on the winding phases
  (node 0 for a solidly earthed phase-to-earth winding; the two nodes of the to side are swapped for clock 6) with
  kVs = u1, u2, XHL and %loadloss from uk and pk, taps from the tap changer and ppm_antifloat=0. An earthing impedance
  is a Reactor from node 4 of the bus (the earthed end of the winding) to earth. The cases use i0 = p0 = 0, since
  OpenDSS places the magnetising branch differently from the pi model of PGM, and only from-side taps (both models
  then keep the ohmic impedance on the to side). The low susceptance that PGM adds to earth on the two phases of a
  phase-to-phase winding (about 1e-8 of sn / uk) is mirrored by a shunt Reactor whose power is added to that side of
  the transformer; without it the transformer powers differ by about 5e-7 relative. Each transformer is first checked
  on its own (open-circuit ratio and short-circuit impedance, so kVs are the winding voltages); node voltages,
  winding powers (without the conductors at earth and at node 4) and the loading (the larger complex sum per side
  divided by sn) are the reference. In a short circuit the winding currents are the reference (with the current of the
  low-susceptance reactor added to its side, as for the powers), and in the OpenDSS netlist each electrode current is
  checked to equal the winding current (KCL; this checks the reference model, not PGM)
"""

import json
import math
import sys
from dataclasses import dataclass, field, replace
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
    # let OpenDSS Kron-reduce the neutral itself (LineCode with Kron=yes) instead of giving it the numpy-reduced matrix
    opendss_kron: bool = False

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


WINDING_CONNECTION = {"a": 1, "b": 2, "c": 3, "ab": 4, "ac": 5, "bc": 6}
BRANCH_SIDE = {"from": 0, "to": 1}


@dataclass
class Transformer1P:
    id: int
    from_node: int
    to_node: int
    from_conn: str  # winding phases: "a", "b", "c" (phase-to-earth) or "ab", "ac", "bc" (phase-to-phase)
    to_conn: str
    u1: float  # winding voltages (V)
    u2: float
    sn: float
    uk: float
    pk: float
    i0: float = 0.0
    p0: float = 0.0
    clock: int = 0
    # earthing impedance of a phase-to-earth winding; None means solidly earthed
    r_g_from: float | None = None
    x_g_from: float | None = None
    r_g_to: float | None = None
    x_g_to: float | None = None
    tap_pos: int | None = None
    tap_nom: int | None = None
    tap_min: int | None = None
    tap_max: int | None = None
    tap_size: float | None = None
    tap_side: str | None = None  # "from" or "to"
    from_status: int = 1
    to_status: int = 1

    def pgm(self):
        result = {
            "id": self.id,
            "from_node": self.from_node,
            "to_node": self.to_node,
            "from_status": self.from_status,
            "to_status": self.to_status,
            "u1": self.u1,
            "u2": self.u2,
            "sn": self.sn,
            "uk": self.uk,
            "pk": self.pk,
            "i0": self.i0,
            "p0": self.p0,
            "from_connection": WINDING_CONNECTION[self.from_conn],
            "to_connection": WINDING_CONNECTION[self.to_conn],
            "clock": self.clock,
        }
        optional = {
            "tap_side": None if self.tap_side is None else BRANCH_SIDE[self.tap_side],
            "tap_pos": self.tap_pos,
            "tap_min": self.tap_min,
            "tap_max": self.tap_max,
            "tap_nom": self.tap_nom,
            "tap_size": self.tap_size,
            "r_grounding_from": self.r_g_from,
            "x_grounding_from": self.x_g_from,
            "r_grounding_to": self.r_g_to,
            "x_grounding_to": self.x_g_to,
        }
        result.update({k: v for k, v in optional.items() if v is not None})
        return result

    def taps(self):
        """OpenDSS taps (per unit of the winding voltage) of the from and to windings, as PGM applies them"""
        if self.tap_pos is None or self.tap_size is None:
            return 1.0, 1.0
        nom = self.tap_nom or 0
        tap_min, tap_max = self.tap_min, self.tap_max
        # the tap direction from the provided limits: reversed when tap_max < tap_min, or with only tap_min above
        # tap_nom, or with only tap_max below tap_nom; otherwise (also without limits) normal
        if tap_min is not None and tap_max is not None:
            reversed_taps = tap_max < tap_min
        elif tap_min is not None:
            reversed_taps = tap_min > nom
        elif tap_max is not None:
            reversed_taps = tap_max < nom
        else:
            reversed_taps = False
        direction = -1 if reversed_taps else 1
        # a provided limit limits the tap position in that direction; a limit not provided leaves its side unbounded
        low, high = (tap_max, tap_min) if reversed_taps else (tap_min, tap_max)
        tap_pos = self.tap_pos
        if high is not None:
            tap_pos = min(tap_pos, high)
        if low is not None:
            tap_pos = max(tap_pos, low)
        offset = direction * (tap_pos - nom) * self.tap_size
        if self.tap_side == "to":
            return 1.0, (self.u2 + offset) / self.u2
        return (self.u1 + offset) / self.u1, 1.0

    def dss_parameters(self):
        """the OpenDSS transformer definition without the buses"""
        r_pct = 100.0 * self.pk / self.sn
        x_pct = 100.0 * math.sqrt(self.uk**2 - (self.pk / self.sn) ** 2)
        imag_pct = 100.0 * math.sqrt(self.i0**2 - (self.p0 / self.sn) ** 2)
        t1, t2 = self.taps()
        return (
            f"conns=[wye wye] kVs=[{self.u1 / 1e3!r} {self.u2 / 1e3!r}] kVAs=[{self.sn / 1e3!r} {self.sn / 1e3!r}]"
            f" XHL={x_pct!r} %loadloss={r_pct!r} %noloadloss={100.0 * self.p0 / self.sn!r} %imag={imag_pct!r}"
            f" ppm_antifloat=0 taps=[{t1!r} {t2!r}]"
        )

    def low_susceptance(self, bus, conn, u_rated):
        """PGM holds the common mode of a phase-to-phase winding with a small admittance to earth, -j b [[1, 1], [1, 1]]
        on its two phases (in siemens, b = 0.5e-8 sn / uk / (u_rated / sqrt(3))^2, the analogue of the zero-sequence
        low susceptance of transformer). The reference includes it as a shunt reactor with the X matrix of
        b ([[1, 1], [1, 1]] + eps I) (eps = 1e-6, as the matrix itself is singular), so it is exact to 1e-6 of this
        already negligible admittance."""
        b = 0.5e-8 * self.sn / self.uk / (u_rated / math.sqrt(3.0)) ** 2
        eps = 1e-6
        x = np.array([[1.0 + eps, -1.0], [-1.0, 1.0 + eps]]) / (b * eps * (2.0 + eps))
        nodes = "".join(f".{'abc'.index(p) + 1}" for p in conn)
        return (
            f"New Reactor.lo{self.id}{bus} phases=2 bus1={bus}{nodes} rmatrix={fmt_matrix(np.zeros((2, 2)))}"
            f" xmatrix={fmt_matrix(x)}"
        )

    def dss_lines(self, u_from, u_to):
        """the OpenDSS transformer, its earth electrodes and the low susceptance of its phase-to-phase windings; u_from
        and u_to are the u_rated of the nodes"""
        buses = []
        shunts = []
        for side, node, conn, r_g, x_g, status, u_rated in (
            ("from", self.from_node, self.from_conn, self.r_g_from, self.x_g_from, self.from_status, u_from),
            ("to", self.to_node, self.to_conn, self.r_g_to, self.x_g_to, self.to_status, u_to),
        ):
            bus = f"n{node}" if status else f"open_t{self.id}_{side}"
            if len(conn) == 2:
                assert r_g is None, "a phase-to-phase winding has no earthed terminal"
                assert x_g is None, "a phase-to-phase winding has no earthed terminal"
                terminals = ["abc".index(p) + 1 for p in conn]
                if status:
                    shunts.append(self.low_susceptance(bus, conn, u_rated))
            elif r_g is None and x_g is None:
                terminals = ["abc".index(conn) + 1, 0]
            else:
                # node 4 of the bus is the earthed end of the winding, connected to earth through the electrode; OpenDSS
                # may treat a Reactor with X = 0 specially, so the cases use x > 0
                assert x_g is not None, "use x_grounding > 0 for an electrode"
                assert x_g > 0.0, "use x_grounding > 0 for an electrode"
                terminals = ["abc".index(conn) + 1, 4]
                shunts.append(
                    f"New Reactor.g{self.id}{side} phases=1 bus1={bus}.4 bus2={bus}.0 R={r_g or 0.0!r} X={x_g or 0.0!r}"
                )
            if side == "to" and self.clock == 6:
                # reversed polarity of the to-side winding
                terminals.reverse()
            buses.append(bus + "".join(f".{t}" for t in terminals))
        return [
            f"New Transformer.t{self.id} phases=1 windings=2 buses=[{buses[0]} {buses[1]}] {self.dss_parameters()}",
            *shunts,
        ]


# OpenDSS replaces r = 0 of a Fault by 1e-4 ohm, so a bolted fault gets this resistance instead
DSS_BOLTED_FAULT_R = 1e-13


@dataclass
class Fault:
    id: int
    node: int
    fault_type: str
    fault_phase: str
    r_f: float
    x_f: float = 0.0

    def dss_r(self):
        return self.r_f if self.r_f > 0.0 else DSS_BOLTED_FAULT_R


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
    # u_rated per node where it differs from u_rated (the source node always has u_rated)
    node_u_rated: dict[int, float] = field(default_factory=dict)
    transformers: list = field(default_factory=list)

    def u_node(self, node):
        return self.node_u_rated.get(node, self.u_rated)


# ---------------------------------------------------------------------------------------------------------------
# PGM input


def pgm_input(case: Case, fault: Fault | None):
    data = {
        "node": [{"id": n, "u_rated": case.u_node(n)} for n in case.nodes],
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
    if case.transformers:
        data["single_phase_transformer"] = [transformer.pgm() for transformer in case.transformers]
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
        if line.opendss_kron:
            assert line.neutral
            assert not line.swer
            z4 = line.z_full() / line.length
            n = z4.shape[0]
            lines.append(
                f"New LineCode.lc{line.id} nphases={n} units=km rmatrix={fmt_matrix(z4.real)}"
                f" xmatrix={fmt_matrix(z4.imag)} cmatrix={fmt_matrix(np.zeros((n, n)))} Neutral={n} Kron=yes"
            )
            lines.append(
                f"New Line.l{line.id} phases={len(line.phases)} bus1={bus1}{nodes} bus2={bus2}{nodes}"
                f" linecode=lc{line.id} length={line.length!r} units=km"
                f" cmatrix={fmt_matrix(c_nf / line.length)}"
            )
            continue
        lines.append(
            f"New Line.l{line.id} phases={len(line.phases)} bus1={bus1}{nodes} bus2={bus2}{nodes} length=1 units=km"
            f" rmatrix={fmt_matrix(z.real)} xmatrix={fmt_matrix(z.imag)} cmatrix={fmt_matrix(c_nf)}"
        )
    electrode_buses = [
        node
        for t in case.transformers
        for node, r_g, x_g in ((t.from_node, t.r_g_from, t.x_g_from), (t.to_node, t.r_g_to, t.x_g_to))
        if r_g is not None or x_g is not None
    ]
    assert len(electrode_buses) == len(set(electrode_buses)), "one electrode (node 4) per bus"
    for transformer in case.transformers:
        lines += transformer.dss_lines(case.u_node(transformer.from_node), case.u_node(transformer.to_node))
    if with_loads:
        for load in case.loads:
            if not load.status:
                continue
            kv_phase = case.u_node(load.node) / math.sqrt(3.0) / 1e3
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
        assert fault.x_f == 0.0, "OpenDSS faults are resistive"
        if fault.fault_type == "single_phase_to_ground":
            k = "abc".index(fault.fault_phase) + 1
            lines.append(f"New Fault.f{fault.id} phases=1 bus1=n{fault.node}.{k} r={fault.dss_r()!r}")
        elif fault.fault_type == "two_phase":
            k1, k2 = ("abc".index(p) + 1 for p in fault.fault_phase)
            lines.append(
                f"New Fault.f{fault.id} phases=1 bus1=n{fault.node}.{k1} bus2=n{fault.node}.{k2} r={fault.dss_r()!r}"
            )
        elif fault.fault_type == "three_phase":
            # three phase-to-earth resistances r_f, as PGM models a three-phase fault with impedance
            lines.append(
                f"New Fault.f{fault.id} phases=3 bus1=n{fault.node}.1.2.3 bus2=n{fault.node}.0.0.0 r={fault.dss_r()!r}"
            )
        else:
            raise NotImplementedError(fault.fault_type)
    bases = ", ".join(repr(u / 1e3) for u in sorted({case.u_node(n) for n in case.nodes}, reverse=True))
    lines += [f"Set VoltageBases=[{bases}]", "CalcVoltageBases", "Set Tolerance=1e-12", "Set MaxIterations=200"]
    if fault is not None:
        lines.append("Set Mode=Snapshot")
    lines.append("Solve")
    return lines


def run_dss(script, case: Case | None = None):
    dss.Basic.ClearAll()
    for command in script:
        dss.Text.Command(command)
    assert dss.Solution.Converged(), "OpenDSS did not converge"
    if case is not None:
        check_opendss_kron_lines(case)


def check_opendss_kron_lines(case: Case):
    """a line with opendss_kron is reduced by OpenDSS itself, to the phases only, and equals the numpy Kron"""
    for line in case.lines:
        if not line.opendss_kron:
            continue
        n = len(line.phases)
        dss.Circuit.SetActiveElement(f"Line.l{line.id}")
        assert dss.CktElement.NumPhases() == n, f"OpenDSS line {line.id} has {dss.CktElement.NumPhases()} phases"
        assert dss.CktElement.NumConductors() == n, f"OpenDSS line {line.id} still has a neutral conductor"
        assert dss.CktElement.NumTerminals() == 2
        assert dss.Lines.Phases() == n
        dss.LineCodes.Name(f"lc{line.id}")
        assert dss.LineCodes.Phases() == n, f"OpenDSS did not Kron-reduce the LineCode of line {line.id}"
        dss.Circuit.SetActiveElement(f"Line.l{line.id}")
        r = np.array(dss.Lines.RMatrix()).reshape(n, n) * line.length
        x = np.array(dss.Lines.XMatrix()).reshape(n, n) * line.length
        expected = line.z_phase()
        dev = np.max(np.abs(r + 1j * x - expected))
        assert dev < 1e-12 * np.max(np.abs(expected)), f"OpenDSS Kron differs from numpy for line {line.id}: {dev}"
        print(f"  OpenDSS-reduced line {line.id}: {n} phases, matches numpy Kron (max dev {dev:.1e} ohm)")


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


def check_transformer(transformer: Transformer1P):
    """sanity checks of the OpenDSS single-phase transformer on its own, fed with the from-side winding voltage:

    - open circuit: the to-side winding voltage is u2 (times the tap ratio), so kVs are the winding voltages
    - short circuit: the impedance seen from the from side is uk * u1^2 / sn (with pk / sn as the resistive part),
      so the per-unit base of a phases=1 transformer is the winding voltage (and not, e.g., the voltage / sqrt(3))
    """
    assert transformer.i0 == 0.0, "the OpenDSS cases use i0 = p0 = 0, see the script"
    assert transformer.p0 == 0.0, "the OpenDSS cases use i0 = p0 = 0, see the script"
    assert transformer.tap_side in (None, "from"), "a to-side tap changes the ohmic impedance differently"
    t1, t2 = transformer.taps()
    u1 = transformer.u1
    params = transformer.dss_parameters()
    results = {}
    for test, secondary in (("open", ""), ("short", "New Fault.short phases=1 bus1=sec.1 r=1e-9")):
        dss.Basic.ClearAll()
        for command in (
            "Clear",
            f"Set DefaultBaseFrequency={FREQ:g}",
            f"New Circuit.t phases=1 basekv={u1 / 1e3!r} pu=1 angle=0 bus1=src.1 Z1=[0, 1e-9] Z0=[0, 1e-9]",
            f"New Transformer.t phases=1 windings=2 buses=[src.1.0 sec.1.0] {params}",
            secondary,
            "Set VoltageBases=[1]",
            "CalcVoltageBases",
            "Set Tolerance=1e-12",
            "Solve",
        ):
            if command:
                dss.Text.Command(command)
        assert dss.Solution.Converged()
        dss.Circuit.SetActiveBus("src")
        v1 = complex(*dss.Bus.Voltages()[:2])
        dss.Circuit.SetActiveBus("sec")
        v2 = complex(*dss.Bus.Voltages()[:2])
        dss.Circuit.SetActiveElement("Transformer.t")
        i1 = complex(*dss.CktElement.Currents()[:2])
        results[test] = (v1, v2, i1)
    v1, v2, _ = results["open"]
    assert abs(abs(v1) / u1 - 1.0) < 1e-9, f"the 1-phase source is not at the winding voltage: {abs(v1)}"
    ratio = v2 / v1
    expected_ratio = transformer.u2 * t2 / (u1 * t1)
    assert abs(ratio / expected_ratio - 1.0) < 1e-9, f"open-circuit ratio {ratio} instead of {expected_ratio}"
    v1, _, i1 = results["short"]
    z = v1 / i1
    r_pu = transformer.pk / transformer.sn
    z_to = complex(r_pu, math.sqrt(transformer.uk**2 - r_pu**2)) * transformer.u2**2 / transformer.sn
    expected_z = z_to * (u1 * t1 / (transformer.u2 * t2)) ** 2
    assert abs(z / expected_z - 1.0) < 1e-6, f"short-circuit impedance {z} instead of {expected_z}"
    print(
        f"  transformer {transformer.id}: open-circuit ratio and short-circuit impedance agree"
        f" (dev {abs(ratio / expected_ratio - 1.0):.1e}, {abs(z / expected_z - 1.0):.1e})"
    )


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
        for branch in case.lines + case.transformers:
            if branch.from_status and branch.to_status:
                for a, b in ((branch.from_node, branch.to_node), (branch.to_node, branch.from_node)):
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
    result = {"node": nodes, "asym_line": lines}
    if case.transformers:
        result["single_phase_transformer"] = [transformer_pf_reference(t) for t in case.transformers]
    return result


def transformer_pf_reference(transformer: Transformer1P):
    """per-phase powers of both windings (the conductors at earth and at the electrode are dropped) and the loading"""
    s_from, s_to = terminal_values(f"Transformer.t{transformer.id}", dss.CktElement.Powers, 2)
    # PGM includes the low susceptance of a phase-to-phase winding in the branch, so its power belongs to that side
    for k, (node, conn) in enumerate(
        ((transformer.from_node, transformer.from_conn), (transformer.to_node, transformer.to_conn))
    ):
        if len(conn) == 2 and (transformer.from_status, transformer.to_status)[k]:
            s_low, _ = terminal_values(f"Reactor.lo{transformer.id}n{node}", dss.CktElement.Powers, 2)
            (s_from if k == 0 else s_to)[:] += s_low
    s_from = s_from * 1e3 if transformer.from_status else np.zeros(3, dtype=complex)
    s_to = s_to * 1e3 if transformer.to_status else np.zeros(3, dtype=complex)
    return {
        "id": transformer.id,
        "p_from": [float(v.real) for v in s_from],
        "q_from": [float(v.imag) for v in s_from],
        "p_to": [float(v.real) for v in s_to],
        "q_to": [float(v.imag) for v in s_to],
        "loading": float(max(abs(np.sum(s_from)), abs(np.sum(s_to))) / transformer.sn),
    }


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
    result = {"node": nodes, "asym_line": lines}
    if case.transformers:
        result["single_phase_transformer"] = [transformer_sc_reference(t, i_f) for t in case.transformers]
    result["fault"] = [{"id": fault.id, "i_f": i_f_mag, "i_f_angle": i_f_angle}]
    return result


def transformer_sc_reference(transformer: Transformer1P, i_f):
    """per-phase winding currents of both sides (the conductors at earth and at the electrode are dropped)"""
    i_from, i_to = terminal_values(f"Transformer.t{transformer.id}", dss.CktElement.Currents, 2)
    for k, (node, conn) in enumerate(
        ((transformer.from_node, transformer.from_conn), (transformer.to_node, transformer.to_conn))
    ):
        if len(conn) == 2 and (transformer.from_status, transformer.to_status)[k]:
            # PGM includes the low susceptance of a phase-to-phase winding in the branch, as for the powers
            i_low, _ = terminal_values(f"Reactor.lo{transformer.id}n{node}", dss.CktElement.Currents, 2)
            (i_from if k == 0 else i_to)[:] += i_low
    check_electrodes(transformer, i_f)
    if not transformer.from_status:
        i_from = np.zeros(3, dtype=complex)
    if not transformer.to_status:
        i_to = np.zeros(3, dtype=complex)
    i_from_mag, i_from_angle = mag_angle(i_from)
    i_to_mag, i_to_angle = mag_angle(i_to)
    return {
        "id": transformer.id,
        "i_from": i_from_mag,
        "i_from_angle": i_from_angle,
        "i_to": i_to_mag,
        "i_to_angle": i_to_angle,
    }


def check_electrodes(transformer: Transformer1P, i_f):
    """a consistency check of the OpenDSS netlist (KCL), not of PGM: the current from node 4 into the electrode equals
    the current into the winding at its phase terminal (the winding has no other path), so the winding current that
    PGM is compared with is also the electrode current"""
    for side, terminal, r_g, x_g, status in (
        ("from", 0, transformer.r_g_from, transformer.x_g_from, transformer.from_status),
        ("to", 1, transformer.r_g_to, transformer.x_g_to, transformer.to_status),
    ):
        if (r_g is None and x_g is None) or not status:
            continue
        dss.Circuit.SetActiveElement(f"Transformer.t{transformer.id}")
        raw = dss.CktElement.Currents()
        order = dss.CktElement.NodeOrder()
        winding = {order[idx]: complex(raw[2 * idx], raw[2 * idx + 1]) for idx in range(2 * terminal, 2 * terminal + 2)}
        phase = "abc".index(transformer.from_conn if side == "from" else transformer.to_conn) + 1
        dss.Circuit.SetActiveElement(f"Reactor.g{transformer.id}{side}")
        i_electrode = complex(*dss.CktElement.Currents()[:2])
        dev = abs(i_electrode - winding[phase])
        assert dev < 1e-9 * abs(i_electrode) + 1e-9, f"electrode of transformer {transformer.id} {side}: {dev}"
        assert abs(winding[4] + winding[phase]) < 1e-9 * abs(i_electrode) + 1e-9
        i_fault = abs(i_f[phase - 1])
        print(
            f"  electrode of transformer {transformer.id} ({side}): {abs(i_electrode):.6g} A, equal to the winding"
            f" current (dev {dev:.1e} A); fault current on phase {'abc'[phase - 1]}: {i_fault:.6g} A"
        )


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
        if case.transformers:
            component_ids["single_phase_transformer"] = [t.id for t in case.transformers]
        write_json(
            target / f"{output_type}.json",
            output_dataset(output_type, {k: [{"id": i} for i in v] for k, v in component_ids.items()}),
        )
        return

    check_kron_against_opendss(case)
    for transformer in case.transformers:
        check_transformer(transformer)

    if case.calculation == "power_flow":
        script = dss_script(case, None, 1.0, with_loads=True)
        run_dss(script, case)
        write_text(target / "test.dss", "\n".join(script) + "\n")
        write_json(target / "asym_output.json", output_dataset("asym_output", pf_reference(case)))
        return

    scripts = []
    references = []
    for fault in case.faults:
        script = dss_script(case, fault, 1.1, with_loads=False)
        run_dss(script, case)
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
# |actual - expected| / (atol + rtol |expected|) is below 0.1 (a 10x margin); see reduced_phase_compare_pgm.py.
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


# u_rated of a node with 240 V phase-to-earth
LV_U = 240.0 * math.sqrt(3.0)


def swer_chain_case():
    """22 kV abc -> isolating transformer ab -> a with an earth electrode -> SWER line -> distribution transformers"""
    return Case(
        path="power_flow/single-phase-transformer/swer-chain",
        u_rated=SWER_U,
        nodes=[1, 2, 3, 4, 5, 6],
        node_u_rated={5: LV_U, 6: LV_U},
        lines=[Line(101, 1, 2, "abc", 2.0), Line(102, 3, 4, "a", 20.0, swer=True)],
        transformers=[
            Transformer1P(151, 2, 3, "ab", "a", 22e3, 12.7e3, 200e3, 0.05, 2000.0, r_g_to=5.0, x_g_to=0.1),
            Transformer1P(
                152,
                4,
                5,
                "a",
                "a",
                12.7e3,
                240.0,
                25e3,
                0.04,
                400.0,
                tap_pos=1,
                tap_nom=0,
                tap_min=-2,
                tap_max=2,
                tap_size=0.025 * 12.7e3,
                tap_side="from",
            ),
            Transformer1P(153, 3, 6, "a", "a", 12.7e3, 240.0, 10e3, 0.035, 180.0),
        ],
        loads=[
            Load(301, 3, "asym_load", "const_pq", [30e3, 0.0, 0.0], [10e3, 0.0, 0.0]),
            Load(302, 4, "asym_load", "const_pq", [40e3, 0.0, 0.0], [12e3, 0.0, 0.0]),
            Load(303, 5, "asym_load", "const_pq", [12e3, 0.0, 0.0], [4e3, 0.0, 0.0]),
            Load(304, 5, "asym_load", "const_y", [6e3, 0.0, 0.0], [2e3, 0.0, 0.0]),
            Load(305, 6, "asym_load", "const_pq", [5e3, 0.0, 0.0], [1.5e3, 0.0, 0.0]),
            Load(306, 6, "asym_load", "const_y", [3e3, 0.0, 0.0], [1e3, 0.0, 0.0]),
        ],
        # limited by the iterative_current convergence and by near-zero powers of phase c of line 101; no tolerance with
        # atol <= 1e-6 passes with a 10x margin (rtol 1e-7, atol 1e-6 gives 0.103)
        params=pf_params(["newton_raphson", "iterative_current"], 1e-7, 1e-5),
    )


def swer_chain_short_circuit_case():
    """the swer-chain network without loads; the earth electrode carries the current of a SWER fault"""
    return replace(
        swer_chain_case(),
        path="short_circuit/single-phase-transformer-swer",
        loads=[],
        calculation="short_circuit",
        params=SC_PARAMS,
        faults=[
            Fault(601, 4, "single_phase_to_ground", "a", 0.0),
            Fault(601, 4, "single_phase_to_ground", "a", 10.0),
            Fault(601, 5, "single_phase_to_ground", "a", 0.0),
            # a three-phase fault in a network with reduced-phase nodes is solved asymmetrically (M4)
            Fault(601, 2, "three_phase", "abc", 0.5),
        ],
    )


def phase_to_phase_case(path, all_const_y, methods):
    """11 kV abc -> ab -> ab (415 V, floating common mode) and bc -> a with clock 6 (240 V)"""
    if all_const_y:
        loads = [
            Load(301, 3, "asym_load", "const_y", [20e3, 12e3, 0.0], [5e3, 3e3, 0.0]),
            Load(302, 4, "asym_load", "const_y", [15e3, 0.0, 0.0], [4e3, 0.0, 0.0]),
        ]
    else:
        # the common mode of n3 is set by the loads only; OpenDSS cannot solve it when they are mostly const_pq (its
        # fixed-point iteration ends on the constant-impedance fallback of the loads, or diverges), so most of the load
        # at n3 is const_y, with a small const_pq part
        loads = [
            Load(301, 3, "asym_load", "const_y", [19e3, 11.5e3, 0.0], [4.75e3, 2.875e3, 0.0]),
            Load(303, 3, "asym_load", "const_pq", [1e3, 0.5e3, 0.0], [0.25e3, 0.125e3, 0.0]),
            Load(302, 4, "asym_load", "const_pq", [15e3, 0.0, 0.0], [4e3, 0.0, 0.0]),
        ]
    return Case(
        path=path,
        u_rated=11e3,
        nodes=[1, 2, 3, 4],
        node_u_rated={3: 415.0, 4: LV_U},
        lines=[Line(101, 1, 2, "abc", 1.0)],
        transformers=[
            Transformer1P(151, 2, 3, "ab", "ab", 11e3, 415.0, 50e3, 0.04, 600.0),
            Transformer1P(152, 2, 4, "bc", "a", 11e3, 240.0, 25e3, 0.04, 350.0, clock=6),
        ],
        loads=loads,
        params=pf_params(methods, 1e-8, 1e-8),
    )


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
        # the only case where OpenDSS Kron-reduces the neutral itself (a LineCode with Kron=yes)
        path="power_flow/asym-line/two-phase-bc-neutral-opendss-kron",
        u_rated=10e3,
        nodes=[1, 2],
        lines=[Line(101, 1, 2, "bc", 4.0, neutral=True, opendss_kron=True)],
        loads=[Load(301, 2, "asym_load", "const_pq", [0.0, 180e3, 90e3], [0.0, 50e3, 20e3])],
        params=pf_params(["newton_raphson", "iterative_current"], 1e-7, 1e-8),
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
    swer_chain_case(),
    swer_chain_short_circuit_case(),
    # the floating common mode of n3 makes iterative_current diverge with const_pq loads, so only newton_raphson
    phase_to_phase_case(
        "power_flow/single-phase-transformer/phase-to-phase", all_const_y=False, methods="newton_raphson"
    ),
    phase_to_phase_case(
        "power_flow/single-phase-transformer/phase-to-phase-linear",
        all_const_y=True,
        methods=["linear", "linear_current"],
    ),
    Case(
        path="power_flow/single-phase-transformer/missing-phase-for-winding",
        u_rated=SWER_U,
        nodes=[1, 2, 3],
        lines=[Line(101, 1, 2, "a", 5.0)],
        transformers=[Transformer1P(151, 2, 3, "ab", "a", 22e3, 12.7e3, 200e3, 0.05, 2000.0)],
        params={
            "calculation_method": "newton_raphson",
            "rtol": 1e-6,
            "atol": 1e-4,
            "experimental_features": "enabled",
            "raises": {
                "raises": "MissingPhaseForWinding",
                "reason": "The ab winding of the single-phase transformer is at a node that only has phase a.",
            },
        },
        raises_only=("asym_output", "MissingPhaseForWinding"),
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
