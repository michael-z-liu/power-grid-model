# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""Generator for the validation cases on two Australian rural SWER networks (SMR8 and KLO14), with OpenDSS as the
reference solver.

The networks are the medium-voltage parts of Network_1_Rural_SMR8.xlsx and Network_2_Rural_KLO14.xlsx of the
Team-Nando MV-LV-Networks repository (BSD 3-Clause), at a pinned commit, as modelled by its backend.py. The xlsx files
are downloaded from the pinned commit (or read from --source-dir) and their sha256 is always checked; they are never
kept in this repository. This script contains no network data.

Run from the repository root (opendssdirect.py and openpyxl are not project dependencies, so they are added on the fly):

    uv run --no-project --with opendssdirect.py --with numpy --with openpyxl \\
        python tests/data/generators/australian_swer_networks.py <root> [--source-dir DIR] [--check-only] [networks]

<root> is the tests/data directory to write to; the networks (smr8, klo14) default to both. --check-only runs the
self-checks and the power flow at minimum and maximum demand, prints a summary and writes nothing.

Modelling choices (derived from backend.py; the deviations are listed):

- IDs: an MV node keeps its AusNet Node_ID; the 66 kV source bus is 1_000_000; the LV busbar of distribution
  transformer E is 2_000_000 + E; a line, transformer, capacitor or the zone substation (Element_ID E) is
  3_000_000 + E; the load of distribution transformer E is 4_000_000 + E; the source is 5_000_000 and a fault
  6_000_000. Every OpenDSS bus is n{node id}; Line.l{id}, Transformer.t{id}, Capacitor.c{id}, Load.ld{id}_{a,b,c}
- source: 66 kV, Z1 = 0.52824 + j2.113 ohm and Z0 = 0.59157 + j1.7 ohm (backend's effective edit). PGM has a real
  z01_ratio, so Z0 gets Z1's angle in both tools (OpenDSS Z0 = z01_ratio Z1); this is irrelevant, because the zone
  substation primary is delta and carries no zero-sequence current
- zone substation (66/22.9 kV SMR8, 66/22.7 kV KLO14, 100 MVA, Dyn): a PGM transformer delta/wye_n clock 1, which is
  OpenDSS's default delta-wye (the LV side lags by 30 degrees)
- voltage regulators (22/22 kV): fixed at the neutral tap as a YNyn0 bank on the rating with backend's xhl and
  %loadloss. backend builds three 12.7/1.27 kV autotransformers under regcontrol starting at tap 0, which at the
  neutral tap are a 1:1 connection plus leakage impedance. KLO14's PLENTY_VLY_RISE (22/22 kV, 200 kVA) is a regulator
  too, as in backend (equal kVs)
- ISO transformers (22/12.7 kV): a single_phase_transformer from the phase pair of Conn_Type (R, W, B = a, b, c) to
  phase a of the SWER network, clock 6 if Conn_Type lists the pair in reverse order; backend's open third winding is
  dropped. The OpenDSS buses are written literally from Conn_Type, so OpenDSS independently checks the clock-6 mapping
- distribution transformers, with the LV network replaced by one lumped load on a new LV busbar:
  - three-phase (22/0.433 kV Dyn, Conn_Type rotation ignored, as backend): PGM transformer delta/wye_n clock 1
  - split-phase (22 kV phase-to-phase, 250-0-250 V): centre-tapped single_phase_transformer with u2 = 500 V and the
    Kersting default split of uk and pk; OpenDSS gets the matching XHL = XHT and XLT instead of backend's defaults
  - SWER (12.7/0.24 kV): single_phase_transformer a -> a
  - KLO14's 44 phase-to-earth units with a 22/0.433 kV rating: 12.7/0.25 kV a -> a (both ratings taken
    phase-to-earth, which keeps the turns ratio of the data)
  - taps: tap_pos = (wdg1_tap - 1) / 0.025 on the from side with tap_size = u1 / 40 and limits -2..2; OpenDSS gets the
    tap explicitly (backend's maxtap/mintap are ignored, as in backend, where the tap stays at wdg1_tap)
- every transformer has i0 = p0 = 0 (%noloadloss = %imag = 0, ppm_antifloat = 0): PGM places the magnetising branch
  as a pi, OpenDSS at winding 1. uk and pk come from the sheet as uk = sqrt(xhl^2 + loadloss^2) / 100 and
  pk = loadloss / 100 sn, and the OpenDSS XHL and %loadloss are computed back from uk and pk
- low susceptance: PGM holds the zero sequence of a delta winding with -j 1e-8 sn / uk (per unit), and the common mode
  of a phase-to-phase single-phase winding likewise; both are mirrored by shunt Reactors (Reactor.lo{id}n{node})
- lines: rows named Delete are dropped. A three-phase line is a PGM line with r1, x1, c1, r0, x0, c0 from the
  LineCode (b in uS/km). A one-phase line is an asym_line on phase a with Z1 and C1 only, because that is how OpenDSS
  solves a 1-phase LineCode (r0, x0, b0 are ignored), so the SWER sections have no earth-return reactance, as in
  backend. All 22 kV and 12.7 kV SWER nodes have u_rated 22 kV (the SWER convention of power-grid-model)
- capacitors: a shunt with b1 = b0 = kvar / kV^2 (grounded wye)
- loads: one const_pq asym_load per distribution transformer at pf 0.95 lagging, f times its kVA: f = 0.10 at minimum
  and MAX_DEMAND at maximum demand; at maximum the SWER loads behind an ISO get the coincidence factor
  min(1, S_iso / (MAX_DEMAND sum S_tx)). Split per phase: three-phase P / 3 each, split-phase P / 2 on each leg,
  phase-to-earth all on phase a. OpenDSS: one 1-phase wye Load per non-zero phase (model=1, kV = u_rated / sqrt(3))
- OpenDSS tolerance 1e-11 (1e-12 stalls at round-off on SMR8) and MaxIterations=200
"""

import argparse
import hashlib
import math
import tempfile
import urllib.request
from collections import Counter, defaultdict, deque
from dataclasses import dataclass, field
from decimal import Decimal
from pathlib import Path

import numpy as np
import opendssdirect as dss
from openpyxl import load_workbook
from reduced_phase_opendss_cases import (
    CENTER_TAPPED,
    FREQ,
    Case,
    Load,
    Transformer1P,
    check_center_tapped_transformer,
    check_transformer,
    dss_load_elements,
    fmt_matrix,
    short_circuit_impedance,
    source_z,
    terminal_values,
    transformer_pf_reference,
)

OMEGA = 2.0 * math.pi * FREQ

SOURCE_REPO = "https://github.com/Team-Nando/MV-LV-Networks"
SOURCE_COMMIT = "8429f8a19f7b95acca2965cc48f2736880e5b31c"
RAW_URL = "https://raw.githubusercontent.com/Team-Nando/MV-LV-Networks/{commit}/{file}"

# demand as a fraction of the transformer kVA, at pf 0.95 lagging
POWER_FACTOR = 0.95
MIN_DEMAND = 0.10
# pending decision D-M10-1: the 0.60 of the approved plan is infeasible (OpenDSS does not converge and the voltage
# collapses to 0.29 pu on SMR8 and 0.53 pu on KLO14, see the README); 0.20 keeps every regulator within its rating
MAX_DEMAND = 0.20
# self-check 5: the lowest node voltage of every class at either demand
MIN_VOLTAGE_PU = 0.8

# the ID scheme (see the module docstring)
SOURCE_NODE = 1_000_000
LV_NODE_OFFSET = 2_000_000
BRANCH_OFFSET = 3_000_000
LOAD_OFFSET = 4_000_000
SOURCE_ID = 5_000_000
FAULT_ID = 6_000_000

U_SOURCE = 66e3
U_MV = 22e3
U_SWER = 12.7e3
SOURCE_Z1 = complex(0.52824, 2.113)
SOURCE_Z0 = complex(0.59157, 1.7)

DSS_TOLERANCE = 1e-11
DSS_MAX_ITERATIONS = 200
DSS_VOLTAGE_BASES = "[66, 22, 0.433, 0.4157]"

# PGM's transformer_low_susceptance_ratio and the regularisation of the mirror Reactor
LOW_SUSCEPTANCE_RATIO = 1e-8
LOW_SUSCEPTANCE_EPS = 1e-6

PHASE_OF_COLOUR = {"R": "a", "W": "b", "B": "c"}
WINDING_WYE_N = 1
WINDING_DELTA = 2
DSS_CONN = {WINDING_WYE_N: "wye", WINDING_DELTA: "delta"}

# the node classes, with the nominal phase-to-earth voltage (V) of the voltage summary and of self-check 5
CLASS_SOURCE = "66 kV source"
CLASS_MV = "22 kV"
CLASS_SWER = "SWER 12.7 kV"
CLASS_LV_3PH = "LV 433 V three-phase"
CLASS_LV_SPLIT = "LV 250-0-250 V split-phase"
CLASS_LV_SWER = "LV 240 V SWER"
CLASS_LV_1PH = "LV 250 V 22 kV phase-to-earth"
NOMINAL_VOLTAGE = {
    CLASS_SOURCE: U_SOURCE / math.sqrt(3.0),
    CLASS_MV: U_MV / math.sqrt(3.0),
    CLASS_SWER: U_SWER,
    CLASS_LV_3PH: 433.0 / math.sqrt(3.0),
    CLASS_LV_SPLIT: 250.0,
    CLASS_LV_SWER: 240.0,
    CLASS_LV_1PH: 250.0,
}

# distribution transformer kinds
KIND_3PH = "three-phase"
KIND_SPLIT = "split-phase"
KIND_SWER = "SWER"
KIND_1PH_22KV = "22 kV phase-to-earth"

D_M10_1_MESSAGE = (
    "maximum demand MAX_DEMAND = {max_demand} of the transformer kVA is infeasible on {network} (decision D-M10-1 of"
    " the M10 plan: the approved 0.60 does not converge in OpenDSS and collapses the voltage; 0.20 is the recommended"
    " maximum)"
)


@dataclass(frozen=True)
class NetworkSpec:
    name: str
    file: str
    sha256: str
    # expected counts (measured during planning)
    lines_three_phase: int
    lines_one_phase: int
    mv_nodes: int
    independent_cycles: int
    cycle_edges: int
    distribution_transformers: dict[str, int]
    isos: int
    regulators: int
    capacitors: int


NETWORKS = (
    NetworkSpec(
        name="smr8",
        file="Network_1_Rural_SMR8.xlsx",
        sha256="53dbce2455ab1014e1fb4776edd1a0e34c4f231fc54c671b0c8a40dc98d0b03f",
        lines_three_phase=2552,
        lines_one_phase=1059,
        mv_nodes=3624,
        independent_cycles=0,
        cycle_edges=0,
        distribution_transformers={KIND_3PH: 253, KIND_SPLIT: 232, KIND_SWER: 219, KIND_1PH_22KV: 0},
        isos=10,
        regulators=2,
        capacitors=2,
    ),
    NetworkSpec(
        name="klo14",
        file="Network_2_Rural_KLO14.xlsx",
        sha256="19d745c6b6aeedf78b2142d5158b4d7761e9e8b0dd35f4f311156a2187cdc46a",
        lines_three_phase=3259,
        lines_one_phase=565,
        mv_nodes=3832,
        # a 64-element loop, a 5-element loop and a parallel pair of connectors, all three-phase
        independent_cycles=3,
        cycle_edges=71,
        distribution_transformers={KIND_3PH: 137, KIND_SPLIT: 358, KIND_SWER: 161, KIND_1PH_22KV: 44},
        isos=7,
        regulators=3,
        capacitors=2,
    ),
)


def scaled(value, factor: int) -> float:
    """a sheet value (kV, kVA) in SI units, exact in decimal (22.9 kV -> 22900.0 V, not 22900.000000000004)"""
    return float(Decimal(repr(value)) * factor)


def tap_position(wdg1_tap) -> int:
    """the tap position of a distribution transformer: wdg1_tap = 1 + 0.025 tap_pos exactly"""
    tap_pos = round((wdg1_tap - 1.0) / 0.025)
    assert abs(1.0 + 0.025 * tap_pos - wdg1_tap) < 1e-12, f"wdg1_tap {wdg1_tap} is not a multiple of 0.025"
    assert -2 <= tap_pos <= 2, f"wdg1_tap {wdg1_tap} is outside 0.95..1.05"
    return tap_pos


def sheet_uk_pk(row, sn: float) -> tuple[float, float]:
    """uk and pk from the sheet's xhl and loadloss (both in percent on the rating)"""
    return math.hypot(row["xhl"], row["loadloss"]) / 100.0, row["loadloss"] / 100.0 * sn


# ---------------------------------------------------------------------------------------------------------------
# source data


def fetch_sources(specs, source_dir: Path | None, download_dir: Path) -> dict[str, Path]:
    """the path of each network's xlsx: from source_dir, or downloaded from the pinned commit; the sha256 is checked"""
    paths = {}
    for spec in specs:
        if source_dir is not None:
            path = source_dir / spec.file
        else:
            url = RAW_URL.format(commit=SOURCE_COMMIT, file=spec.file)
            print(f"downloading {url}")
            path = download_dir / spec.file
            with urllib.request.urlopen(url, timeout=300) as response:
                path.write_bytes(response.read())
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != spec.sha256:
            raise SystemExit(f"{path}: sha256 {digest} instead of {spec.sha256}")
        print(f"{spec.file}: sha256 {digest} (as pinned)")
        paths[spec.name] = path
    return paths


def read_workbook(path: Path) -> dict[str, list[dict]]:
    """every sheet as a list of rows (dicts by the header row)"""
    workbook = load_workbook(path, read_only=True, data_only=True)
    sheets = {}
    for sheet in workbook.worksheets:
        rows = sheet.iter_rows(values_only=True)
        header = next(rows)
        sheets[sheet.title] = [dict(zip(header, row)) for row in rows]
    workbook.close()
    return sheets


# ---------------------------------------------------------------------------------------------------------------
# network model


@dataclass
class ThreePhaseTransformer:
    """a PGM transformer (zone substation, regulator or three-phase distribution transformer) with only from-side taps
    and i0 = p0 = 0"""

    id: int
    from_node: int
    to_node: int
    u1: float
    u2: float
    sn: float
    uk: float
    pk: float
    winding_from: int
    winding_to: int
    clock: int
    u_rated_from: float
    u_rated_to: float
    tap_pos: int = 0
    tap_min: int = 0
    tap_max: int = 0
    tap_size: float = 0.0

    def pgm(self):
        return {
            "id": self.id,
            "from_node": self.from_node,
            "to_node": self.to_node,
            "from_status": 1,
            "to_status": 1,
            "u1": self.u1,
            "u2": self.u2,
            "sn": self.sn,
            "uk": self.uk,
            "pk": self.pk,
            "i0": 0.0,
            "p0": 0.0,
            "winding_from": self.winding_from,
            "winding_to": self.winding_to,
            "clock": self.clock,
            "tap_side": 0,
            "tap_pos": self.tap_pos,
            "tap_min": self.tap_min,
            "tap_max": self.tap_max,
            "tap_nom": 0,
            "tap_size": self.tap_size,
        }

    def tap_ratio(self) -> float:
        """the OpenDSS tap of the from winding, as PGM applies tap_pos (the direction from the limits)"""
        assert min(self.tap_min, self.tap_max) <= self.tap_pos <= max(self.tap_min, self.tap_max)
        direction = 1 if self.tap_max > self.tap_min else -1
        return (self.u1 + direction * self.tap_pos * self.tap_size) / self.u1

    def low_susceptance(self) -> float | None:
        """the susceptance (S, positive; the admittance is -j b) of PGM's zero-sequence low susceptance on the from
        side, which PGM adds on a side without a zero-sequence path (here: a delta from side)"""
        assert self.winding_to == WINDING_WYE_N, "only a wye_n to side is used"
        if self.winding_from != WINDING_DELTA:
            return None
        return LOW_SUSCEPTANCE_RATIO * self.sn / (self.uk * self.u_rated_from**2)

    def dss_parameters(self) -> str:
        r_pct = 100.0 * self.pk / self.sn
        x_pct = 100.0 * math.sqrt(self.uk**2 - (self.pk / self.sn) ** 2)
        conns = f"{DSS_CONN[self.winding_from]} {DSS_CONN[self.winding_to]}"
        return (
            f"conns=[{conns}] kVs=[{self.u1 / 1e3!r} {self.u2 / 1e3!r}] kVAs=[{self.sn / 1e3!r} {self.sn / 1e3!r}]"
            f" XHL={x_pct!r} %loadloss={r_pct!r} %noloadloss=0 %imag=0 ppm_antifloat=0"
            f" taps=[{self.tap_ratio()!r} 1.0]"
        )

    def mirror_reactor(self, bus: str) -> str | None:
        """the shunt Reactor that mirrors the low susceptance: Y = -j (b / 3) (J + eps I), whose zero-sequence part is
        PGM's -j b and whose positive and negative sequence (-j b eps / 3) is negligible"""
        b = self.low_susceptance()
        if b is None:
            return None
        eps = LOW_SUSCEPTANCE_EPS
        x = 3.0 / (b * eps) * (np.eye(3) - np.ones((3, 3)) / (3.0 + eps))
        return (
            f"New Reactor.lo{self.id}{bus} phases=3 bus1={bus}.1.2.3 rmatrix={fmt_matrix(np.zeros((3, 3)))}"
            f" xmatrix={fmt_matrix(x)}"
        )

    def dss_lines(self) -> list[str]:
        lines = [
            (
                f"New Transformer.t{self.id} phases=3 windings=2 buses=[n{self.from_node} n{self.to_node}]"
                f" {self.dss_parameters()}"
            )
        ]
        mirror = self.mirror_reactor(f"n{self.from_node}")
        if mirror is not None:
            lines.append(mirror)
        return lines

    def primitive_admittance(self) -> np.ndarray:
        """the 6x6 primitive admittance (siemens; from a, b, c, to a, b, c) of the PGM model without the low
        susceptance, written independently of PGM after transformer.hpp: per block (ff, ft, tf, tt)
        Y_abc = A diag(y0, y1, y2) A^-1 with the per-unit sequence admittances of Edge::calc_param_y_sym, scaled with
        1e6 / (u_i u_j) on the node u_rated"""
        u1 = self.u1 * self.tap_ratio()
        k = (u1 / self.u2) / (self.u_rated_from / self.u_rated_to)
        z_ohm = short_circuit_impedance(self.uk, self.pk, self.sn) * self.u2**2 / self.sn
        y = 1.0 / z_ohm * self.u_rated_to**2 / 1e6

        def blocks(ratio: complex):
            return np.array([y / abs(ratio) ** 2, -y / ratio.conjugate(), -y / ratio, y])

        shift = self.clock * math.pi / 6.0
        y1 = blocks(k * np.exp(1j * shift))
        y2 = blocks(k * np.exp(-1j * shift))
        if self.winding_from == WINDING_WYE_N and self.winding_to == WINDING_WYE_N:
            shift_0 = math.pi if self.clock in (2, 6, 10) else 0.0
            y0 = blocks(k * np.exp(1j * shift_0))
        elif self.winding_from == WINDING_DELTA and self.winding_to == WINDING_WYE_N:
            y0 = np.array([0.0, 0.0, 0.0, y])
        else:
            raise NotImplementedError((self.winding_from, self.winding_to))
        a = np.exp(2j * math.pi / 3.0)
        sym = np.array([[1.0, 1.0, 1.0], [1.0, a * a, a], [1.0, a, a * a]])
        sym_inv = np.linalg.inv(sym)
        u = (self.u_rated_from, self.u_rated_to)
        result = np.zeros((6, 6), dtype=complex)
        for block, (i, j) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
            y_abc = sym @ np.diag([y0[block], y1[block], y2[block]]) @ sym_inv
            result[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] = y_abc * 1e6 / (u[i] * u[j])
        return result


@dataclass
class LineSection:
    id: int
    from_node: int
    to_node: int
    phases: int  # 3 or 1
    length: float  # km
    code: str  # the LineCode key, e.g. "46-3ph"


@dataclass
class Network:
    spec: NetworkSpec
    linecodes: dict[str, dict]  # the used LineCodes by key
    lines: list[LineSection]
    node_class: dict[int, str]
    node_u_rated: dict[int, float]
    zone: ThreePhaseTransformer
    regulators: list[ThreePhaseTransformer]
    isos: list[Transformer1P]
    iso_names: dict[int, str]
    iso_buses: dict[int, str]  # the literal OpenDSS primary bus of each ISO
    regulator_names: dict[int, str]
    distribution: list[ThreePhaseTransformer | Transformer1P]
    distribution_kind: dict[int, str]
    distribution_buses: dict[int, list[str]]  # the literal OpenDSS buses of each single-phase distribution transformer
    distribution_kva: dict[int, float]
    capacitors: list[dict]
    iso_of_load: dict[int, int]  # load id -> the ISO that feeds it (SWER loads only)
    iso_factor: dict[int, float]
    iso_swer_kva: dict[int, float]
    loads_min: list[Load] = field(default_factory=list)
    loads_max: list[Load] = field(default_factory=list)

    def case(self) -> Case:
        """the existing generator's Case as an adapter for source_z and dss_load_elements"""
        z1_abs = abs(SOURCE_Z1)
        return Case(
            path=f"australian-swer/{self.spec.name}",
            u_rated=U_SOURCE,
            nodes=sorted(self.node_u_rated),
            lines=[],
            sk=U_SOURCE**2 / z1_abs,
            rx_ratio=SOURCE_Z1.real / SOURCE_Z1.imag,
            z01_ratio=abs(SOURCE_Z0) / z1_abs,
            node_u_rated=self.node_u_rated,
        )

    def three_phase_transformers(self) -> list[ThreePhaseTransformer]:
        return [
            self.zone,
            *self.regulators,
            *(t for t in self.distribution if isinstance(t, ThreePhaseTransformer)),
        ]

    def single_phase_transformers(self) -> list[Transformer1P]:
        return [*self.isos, *(t for t in self.distribution if isinstance(t, Transformer1P))]


def sorted_pair(conn_type: str) -> tuple[str, int]:
    """the PGM phase pair of a two-phase Conn_Type and the clock: 6 if Conn_Type lists the pair in reverse order"""
    phases = [PHASE_OF_COLOUR[c] for c in conn_type]
    pair = "".join(sorted(phases))
    return pair, 0 if "".join(phases) == pair else 6


def dss_phase_nodes(conn_type: str) -> str:
    """the OpenDSS nodes of Conn_Type, in its literal order (R, W, B = 1, 2, 3), as backend's char_to_num"""
    return "".join(f".{'RWB'.index(c) + 1}" for c in conn_type)


def check_unique_element_ids(sheets):
    ids = [
        row["Element_ID"]
        for name in ("lines", "mvtx", "lvtx", "mvcaps", "mv_net_txs")
        for row in sheets[name]
        if name != "lines" or str(row["Element_Name"]).lower() != "delete"
    ]
    assert len(ids) == len(set(ids)), "Element_IDs are not unique"
    assert all(isinstance(i, int) and 0 < i < 1_000_000 for i in ids)


def used_linecodes(sheets, lines: list[LineSection]) -> dict[str, dict]:
    by_key = defaultdict(list)
    for row in sheets["linecodes"]:
        by_key[row["Linecode_ID"]].append(row)
    result = {}
    for key in sorted({line.code for line in lines}):
        rows = by_key[key]
        assert len(rows) == 1, f"LineCode {key} is defined {len(rows)} times"
        row = rows[0]
        assert row["Units"] == "km"
        assert row["Phases"] == int(key.split("-")[1][0])
        result[key] = row
    return result


def parse_lines(sheets) -> list[LineSection]:
    lines = []
    for row in sheets["lines"]:
        if str(row["Element_Name"]).lower() == "delete":
            continue
        phases = row["Phases"]
        expected_phase = "1.2.3" if phases == 3 else 1
        assert phases in (1, 3)
        assert row["Start_Node_Phase"] == expected_phase, row
        assert row["End_Node_Phase"] == expected_phase, row
        assert row["Units"] == "km"
        assert row["Start_Node"] != row["End_Node"]
        assert row["Length"] > 0.0
        lines.append(
            LineSection(
                id=BRANCH_OFFSET + row["Element_ID"],
                from_node=row["Start_Node"],
                to_node=row["End_Node"],
                phases=phases,
                length=row["Length"],
                code=f"{row['Linecode']}-{phases}ph",
            )
        )
    return sorted(lines, key=lambda line: line.id)


def cycle_edges(adjacency: dict[int, list[tuple[int, int]]], root: int) -> set[int]:
    """the edges on a cycle (the non-bridges) of a connected multigraph, by an iterative Tarjan bridge search"""
    order: dict[int, int] = {root: 0}
    low: dict[int, int] = {root: 0}
    bridges = set()
    stack = [(root, -1, iter(adjacency[root]))]
    while stack:
        node, parent_edge, neighbours = stack[-1]
        advanced = False
        for edge, other in neighbours:
            if edge == parent_edge:
                continue
            if other in order:
                low[node] = min(low[node], order[other])
            else:
                order[other] = low[other] = len(order)
                stack.append((other, edge, iter(adjacency[other])))
                advanced = True
                break
        if advanced:
            continue
        stack.pop()
        if stack:
            parent = stack[-1][0]
            low[parent] = min(low[parent], low[node])
            if low[node] > order[parent]:
                bridges.add(parent_edge)
    all_edges = {edge for edges in adjacency.values() for edge, _ in edges}
    return all_edges - bridges


def build_network(spec: NetworkSpec, sheets) -> Network:
    """the section 2 mapping of the plan, with the topology checks (self-check 7)"""
    check_unique_element_ids(sheets)
    lines = parse_lines(sheets)
    linecodes = used_linecodes(sheets, lines)
    n_three = sum(1 for line in lines if line.phases == 3)
    n_one = len(lines) - n_three
    assert (n_three, n_one) == (spec.lines_three_phase, spec.lines_one_phase), (n_three, n_one)

    (zone_row,) = sheets["mv_net_txs"]
    assert zone_row["Bus1"] == "SourceBus"
    assert (zone_row["Connection_Primary"], zone_row["Connection_Secondary"]) == ("delta", "wye")
    assert zone_row["kvs_primary"] == 66
    assert zone_row["kvas_primary"] == zone_row["kvas_secondary"]
    zone_sn = scaled(zone_row["kvas_primary"], 1000)
    zone_uk, zone_pk = sheet_uk_pk(zone_row, zone_sn)
    zone = ThreePhaseTransformer(
        id=BRANCH_OFFSET + zone_row["Element_ID"],
        from_node=SOURCE_NODE,
        to_node=zone_row["Bus2"],
        u1=U_SOURCE,
        u2=scaled(zone_row["kvs_secondary"], 1000),
        sn=zone_sn,
        uk=zone_uk,
        pk=zone_pk,
        winding_from=WINDING_DELTA,
        winding_to=WINDING_WYE_N,
        clock=1,
        u_rated_from=U_SOURCE,
        u_rated_to=U_MV,
    )

    regulators = []
    regulator_names = {}
    isos = []
    iso_names = {}
    iso_buses = {}
    for row in sheets["mvtx"]:
        tid = BRANCH_OFFSET + row["Element_ID"]
        assert row["kvas_primary"] == row["kvas_secondary"]
        sn = scaled(row["kvas_primary"], 1000)
        uk, pk = sheet_uk_pk(row, sn)
        assert row["kvs_primary"] == 22
        if row["kvs_secondary"] == row["kvs_primary"]:
            assert len(row["Conn_Type"]) == 3
            regulators.append(
                ThreePhaseTransformer(
                    id=tid,
                    from_node=row["Bus1"],
                    to_node=row["Bus2"],
                    u1=U_MV,
                    u2=U_MV,
                    sn=sn,
                    uk=uk,
                    pk=pk,
                    winding_from=WINDING_WYE_N,
                    winding_to=WINDING_WYE_N,
                    clock=0,
                    u_rated_from=U_MV,
                    u_rated_to=U_MV,
                )
            )
            regulator_names[tid] = row["Substation_ID"]
            continue
        assert row["kvs_secondary"] == 12.7
        assert len(row["Conn_Type"]) == 2
        pair, clock = sorted_pair(row["Conn_Type"])
        isos.append(Transformer1P(tid, row["Bus1"], row["Bus2"], pair, "a", U_MV, U_SWER, sn, uk, pk, clock=clock))
        iso_names[tid] = row["Substation_ID"]
        iso_buses[tid] = f"n{row['Bus1']}{dss_phase_nodes(row['Conn_Type'])}"
    assert len(isos) == spec.isos
    assert len(regulators) == spec.regulators

    capacitors = []
    for row in sheets["mvcaps"]:
        assert row["phases"] == 3
        assert row["kvs"] == 22
        capacitors.append(
            {
                "id": BRANCH_OFFSET + row["Element_ID"],
                "node": row["Bus1"],
                "kvar": row["kvar"],
                "kv": row["kvs"],
                "b": scaled(row["kvar"], 1000) / scaled(row["kvs"], 1000) ** 2,
            }
        )
    assert len(capacitors) == spec.capacitors

    # graph: (edge id, node a, node b, kind)
    edges = [(line.id, line.from_node, line.to_node, f"line{line.phases}") for line in lines]
    edges.append((zone.id, zone.from_node, zone.to_node, "zone"))
    edges += [(t.id, t.from_node, t.to_node, "regulator") for t in regulators]
    edges += [(t.id, t.from_node, t.to_node, "iso") for t in isos]
    lvtx_rows = sorted(sheets["lvtx"], key=lambda row: row["Element_ID"])
    edges += [
        (BRANCH_OFFSET + row["Element_ID"], row["Bus1"], LV_NODE_OFFSET + row["Element_ID"], "lvtx")
        for row in lvtx_rows
    ]
    mv_nodes = {node for _, a, b, kind in edges for node in (a, b) if kind != "zone" and node < LV_NODE_OFFSET}
    mv_nodes.add(zone.to_node)
    mv_nodes |= {c["node"] for c in capacitors}
    assert len(mv_nodes) == spec.mv_nodes, f"{len(mv_nodes)} MV nodes instead of {spec.mv_nodes}"
    assert all(0 < node < SOURCE_NODE for node in mv_nodes)

    adjacency: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for edge, a, b, _ in edges:
        adjacency[a].append((edge, b))
        adjacency[b].append((edge, a))
    reached = {SOURCE_NODE}
    queue = deque([SOURCE_NODE])
    while queue:
        node = queue.popleft()
        for _, other in adjacency[node]:
            if other not in reached:
                reached.add(other)
                queue.append(other)
    all_nodes = mv_nodes | {SOURCE_NODE} | {LV_NODE_OFFSET + row["Element_ID"] for row in lvtx_rows}
    assert reached == all_nodes, f"{len(all_nodes - reached)} nodes are not reachable from the source"
    n_cycles = len(edges) - len(all_nodes) + 1
    assert n_cycles == spec.independent_cycles, f"{n_cycles} independent cycles"
    on_cycle = cycle_edges(adjacency, SOURCE_NODE)
    kind_of_edge = {edge: kind for edge, _, _, kind in edges}
    assert len(on_cycle) == spec.cycle_edges, f"{len(on_cycle)} edges on a cycle"
    assert all(kind_of_edge[edge] == "line3" for edge in on_cycle), "a cycle has an element other than a 3-phase line"

    # the SWER network behind each ISO: the nodes reached from its secondary through one-phase lines only
    one_phase_adjacency = defaultdict(list)
    for line in lines:
        if line.phases == 1:
            one_phase_adjacency[line.from_node].append(line.to_node)
            one_phase_adjacency[line.to_node].append(line.from_node)
    swer_of_node: dict[int, int] = {}
    for iso in isos:
        assert iso.to_node not in swer_of_node, "two ISOs feed the same SWER network"
        swer_of_node[iso.to_node] = iso.id
        queue = deque([iso.to_node])
        while queue:
            node = queue.popleft()
            for other in one_phase_adjacency[node]:
                if other not in swer_of_node:
                    swer_of_node[other] = iso.id
                    queue.append(other)
                assert swer_of_node[other] == iso.id, "two ISOs feed the same SWER network"
    for line in lines:
        ends_on_swer = (line.from_node in swer_of_node, line.to_node in swer_of_node)
        assert ends_on_swer == ((True, True) if line.phases == 1 else (False, False)), (
            f"line {line.id} mixes the networks"
        )
    assert all(iso.from_node not in swer_of_node for iso in isos)
    assert all(t.from_node not in swer_of_node and t.to_node not in swer_of_node for t in [zone, *regulators])
    assert all(c["node"] not in swer_of_node for c in capacitors), "a capacitor on a SWER node"

    node_class = {SOURCE_NODE: CLASS_SOURCE}
    node_class.update({node: CLASS_SWER if node in swer_of_node else CLASS_MV for node in mv_nodes})
    node_u_rated = {node: U_MV for node in mv_nodes}
    node_u_rated[SOURCE_NODE] = U_SOURCE

    distribution: list[ThreePhaseTransformer | Transformer1P] = []
    distribution_kind = {}
    distribution_buses = {}
    distribution_kva = {}
    for row in lvtx_rows:
        element = row["Element_ID"]
        tid = BRANCH_OFFSET + element
        lv = LV_NODE_OFFSET + element
        bus1 = row["Bus1"]
        conn = row["Conn_Type"]
        assert row["Bus2"] == 0
        assert (row["Connection_Primary"], row["Connection_Secondary"]) == ("delta", "wye")
        assert row["kvas_primary"] == row["kvas_secondary"]
        sn = scaled(row["kvas_primary"], 1000)
        uk, pk = sheet_uk_pk(row, sn)
        tap_pos = tap_position(row["wdg1_tap"])
        kv = (row["kvs_primary"], row["kvs_secondary"])
        on_swer = bus1 in swer_of_node
        distribution_kva[tid] = float(row["kvas_primary"])
        if len(conn) == 3:
            assert kv == (22, 0.433)
            assert not on_swer
            kind = KIND_3PH
            node_class[lv], node_u_rated[lv] = CLASS_LV_3PH, 433.0
            distribution.append(
                ThreePhaseTransformer(
                    id=tid,
                    from_node=bus1,
                    to_node=lv,
                    u1=U_MV,
                    u2=433.0,
                    sn=sn,
                    uk=uk,
                    pk=pk,
                    winding_from=WINDING_DELTA,
                    winding_to=WINDING_WYE_N,
                    clock=1,
                    u_rated_from=U_MV,
                    u_rated_to=433.0,
                    tap_pos=tap_pos,
                    tap_min=-2,
                    tap_max=2,
                    tap_size=U_MV / 40.0,
                )
            )
        else:
            taps = {"tap_side": "from", "tap_pos": tap_pos, "tap_nom": 0, "tap_min": -2, "tap_max": 2}
            if len(conn) == 2:
                assert kv == (22, 0.433)
                assert not on_swer
                kind = KIND_SPLIT
                node_class[lv], node_u_rated[lv] = CLASS_LV_SPLIT, 250.0 * math.sqrt(3.0)
                pair, clock = sorted_pair(conn)
                transformer = Transformer1P(
                    tid,
                    bus1,
                    lv,
                    pair,
                    CENTER_TAPPED,
                    U_MV,
                    500.0,
                    sn,
                    uk,
                    pk,
                    clock=clock,
                    tap_size=U_MV / 40.0,
                    **taps,
                )
                buses = [f"n{bus1}{dss_phase_nodes(conn)}", f"n{lv}.1.0", f"n{lv}.0.2"]
            else:
                assert conn == "R", "a phase-to-earth unit is wired from phase R"
                if kv == (12.7, 0.24):
                    assert on_swer
                    kind = KIND_SWER
                    node_class[lv], u2 = CLASS_LV_SWER, 240.0
                else:
                    assert kv == (22, 0.433)
                    assert not on_swer
                    # decision D-M10-3: both data ratings taken phase-to-earth (the data's turns ratio)
                    kind = KIND_1PH_22KV
                    node_class[lv], u2 = CLASS_LV_1PH, 250.0
                node_u_rated[lv] = u2 * math.sqrt(3.0)
                transformer = Transformer1P(
                    tid, bus1, lv, "a", "a", U_SWER, u2, sn, uk, pk, tap_size=U_SWER / 40.0, **taps
                )
                buses = [f"n{bus1}.1.0", f"n{lv}.1.0"]
            distribution.append(transformer)
            distribution_buses[tid] = buses
        distribution_kind[tid] = kind
    counts = Counter(distribution_kind.values())
    assert {kind: counts[kind] for kind in spec.distribution_transformers} == spec.distribution_transformers, counts

    # demand (decisions D-M10-1 and D-M10-2)
    iso_swer_kva = {iso.id: 0.0 for iso in isos}
    iso_of_load = {}
    for t in distribution:
        if distribution_kind[t.id] == KIND_SWER:
            iso_id = swer_of_node[t.from_node]
            iso_swer_kva[iso_id] += distribution_kva[t.id]
            iso_of_load[t.id - BRANCH_OFFSET + LOAD_OFFSET] = iso_id
    iso_factor = {
        iso.id: min(1.0, iso.sn / 1e3 / (MAX_DEMAND * iso_swer_kva[iso.id])) if iso_swer_kva[iso.id] else 1.0
        for iso in isos
    }

    network = Network(
        spec=spec,
        linecodes=linecodes,
        lines=lines,
        node_class=node_class,
        node_u_rated=node_u_rated,
        zone=zone,
        regulators=regulators,
        isos=isos,
        iso_names=iso_names,
        iso_buses=iso_buses,
        regulator_names=regulator_names,
        distribution=distribution,
        distribution_kind=distribution_kind,
        distribution_buses=distribution_buses,
        distribution_kva=distribution_kva,
        capacitors=capacitors,
        iso_of_load=iso_of_load,
        iso_factor=iso_factor,
        iso_swer_kva=iso_swer_kva,
    )
    network.loads_min = demand(network, MIN_DEMAND, with_iso_factor=False)
    network.loads_max = demand(network, MAX_DEMAND, with_iso_factor=True)
    check_unique_pgm_ids(network)
    print(
        f"{spec.name}: topology as expected: {len(lines)} lines ({n_three} three-phase, {n_one} one-phase),"
        f" {len(mv_nodes)} MV nodes, all reachable, {n_cycles} independent cycles ({len(on_cycle)} edges, all"
        f" three-phase lines), {len(isos)} ISOs each feeding its own SWER network, {len(regulators)} regulators,"
        f" distribution transformers {dict(sorted(counts.items()))}"
    )
    return network


def demand(network: Network, fraction: float, with_iso_factor: bool) -> list[Load]:
    """one const_pq asym_load per distribution transformer: fraction times its kVA at pf 0.95 lagging, times the ISO
    factor for a SWER load at maximum demand"""
    q_per_p = math.sqrt(1.0 - POWER_FACTOR**2) / POWER_FACTOR
    loads = []
    for t in network.distribution:
        load_id = t.id - BRANCH_OFFSET + LOAD_OFFSET
        factor = fraction
        if with_iso_factor and load_id in network.iso_of_load:
            factor *= network.iso_factor[network.iso_of_load[load_id]]
        p = factor * network.distribution_kva[t.id] * 1e3 * POWER_FACTOR
        q = p * q_per_p
        kind = network.distribution_kind[t.id]
        if kind == KIND_3PH:
            p_phase, q_phase = [p / 3.0] * 3, [q / 3.0] * 3
        elif kind == KIND_SPLIT:
            p_phase, q_phase = [p / 2.0, p / 2.0, 0.0], [q / 2.0, q / 2.0, 0.0]
        else:
            p_phase, q_phase = [p, 0.0, 0.0], [q, 0.0, 0.0]
        loads.append(Load(load_id, t.to_node, "asym_load", "const_pq", p_phase, q_phase))
    return loads


def check_unique_pgm_ids(network: Network):
    ids = [
        *network.node_u_rated,
        *(line.id for line in network.lines),
        *(t.id for t in network.three_phase_transformers()),
        *(t.id for t in network.single_phase_transformers()),
        *(c["id"] for c in network.capacitors),
        *(load.id for load in network.loads_min),
        SOURCE_ID,
        FAULT_ID,
    ]
    assert len(ids) == len(set(ids)), "the PGM IDs are not unique"


# ---------------------------------------------------------------------------------------------------------------
# PGM input rows


def pgm_rows(network: Network, loads: list[Load]) -> dict[str, list[dict]]:
    """the PGM input of the network with the given loads, per component, sorted by id"""
    case = network.case()
    lines = []
    asym_lines = []
    for line in network.lines:
        code = network.linecodes[line.code]
        length = line.length
        row: dict[str, int | float] = {
            "id": line.id,
            "from_node": line.from_node,
            "to_node": line.to_node,
            "from_status": 1,
            "to_status": 1,
        }
        i_n = float(min(code["Ampacity1"], code["Ampacity2"]))
        if line.phases == 3:
            row.update(
                {
                    "r1": code["r1"] * length,
                    "x1": code["x1"] * length,
                    "c1": code["b1"] * 1e-6 / OMEGA * length,
                    "tan1": 0.0,
                    "r0": code["r0"] * length,
                    "x0": code["x0"] * length,
                    "c0": code["b0"] * 1e-6 / OMEGA * length,
                    "tan0": 0.0,
                    "i_n": i_n,
                }
            )
            lines.append(row)
        else:
            # OpenDSS solves a 1-phase LineCode with Z1 and C1 only
            row.update(
                {
                    "phases": 1,
                    "r_aa": code["r1"] * length,
                    "x_aa": code["x1"] * length,
                    "c_aa": code["b1"] * 1e-6 / OMEGA * length,
                    "i_n": i_n,
                }
            )
            asym_lines.append(row)
    return {
        "node": [{"id": node, "u_rated": network.node_u_rated[node]} for node in sorted(network.node_u_rated)],
        "line": lines,
        "asym_line": asym_lines,
        "transformer": sorted((t.pgm() for t in network.three_phase_transformers()), key=lambda row: row["id"]),
        "single_phase_transformer": sorted(
            (t.pgm() for t in network.single_phase_transformers()), key=lambda row: row["id"]
        ),
        "shunt": [
            {"id": c["id"], "node": c["node"], "status": 1, "g1": 0.0, "b1": c["b"], "g0": 0.0, "b0": c["b"]}
            for c in sorted(network.capacitors, key=lambda c: c["id"])
        ],
        "source": [
            {
                "id": SOURCE_ID,
                "node": SOURCE_NODE,
                "status": 1,
                "u_ref": 1.0,
                "sk": case.sk,
                "rx_ratio": case.rx_ratio,
                "z01_ratio": case.z01_ratio,
            }
        ],
        "asym_load": [load.pgm() for load in sorted(loads, key=lambda load: load.id)],
    }


# ---------------------------------------------------------------------------------------------------------------
# OpenDSS model


def dss_network(network: Network, voltage_pu: float = 1.0) -> list[str]:
    """the OpenDSS network without loads: the circuit (source), LineCodes, Lines, Capacitors and the transformers with
    their low-susceptance Reactors; voltage_pu = c scales the source voltage and impedance (IEC 60909)"""
    case = network.case()
    z1, z0 = source_z(case)
    if voltage_pu != 1.0:
        z1 *= voltage_pu
        z0 *= voltage_pu
    commands = [
        "Clear",
        f"Set DefaultBaseFrequency={FREQ:g}",
        (
            f"New Circuit.c phases=3 basekv={U_SOURCE / 1e3!r} pu={voltage_pu!r} angle=0 bus1=n{SOURCE_NODE}"
            f" Z1=[{z1.real!r}, {z1.imag!r}] Z0=[{z0.real!r}, {z0.imag!r}]"
        ),
    ]
    for key, code in network.linecodes.items():
        commands.append(
            f"New LineCode.lc{key} nphases={code['Phases']} r1={code['r1']!r} x1={code['x1']!r} b1={code['b1']!r}"
            f" r0={code['r0']!r} x0={code['x0']!r} b0={code['b0']!r} units=km"
        )
    for line in network.lines:
        nodes = ".1.2.3" if line.phases == 3 else ".1"
        commands.append(
            f"New Line.l{line.id} phases={line.phases} bus1=n{line.from_node}{nodes} bus2=n{line.to_node}{nodes}"
            f" linecode=lc{line.code} length={line.length!r} units=km"
        )
    commands += [
        f"New Capacitor.c{c['id']} bus1=n{c['node']}.1.2.3 phases=3 kvar={c['kvar']!r} kV={c['kv']!r}"
        for c in sorted(network.capacitors, key=lambda c: c["id"])
    ]
    transformers: list[ThreePhaseTransformer | Transformer1P] = [
        *network.three_phase_transformers(),
        *network.single_phase_transformers(),
    ]
    for t in sorted(transformers, key=lambda t: t.id):
        if isinstance(t, ThreePhaseTransformer):
            commands += t.dss_lines()
            continue
        if t.id in network.iso_buses:
            buses = [network.iso_buses[t.id], f"n{t.to_node}.1.0"]
        else:
            buses = network.distribution_buses[t.id]
        commands.append(
            f"New Transformer.t{t.id} phases=1 windings={len(buses)} buses=[{' '.join(buses)}] {t.dss_parameters()}"
        )
        if len(t.from_conn) == 2:
            commands.append(t.low_susceptance(f"n{t.from_node}", t.from_conn, network.node_u_rated[t.from_node]))
    return commands


def dss_loads(network: Network, loads: list[Load]) -> list[str]:
    case = network.case()
    return [f"New Load.{name} {definition}" for load in loads for name, definition, _ in dss_load_elements(case, load)]


def dss_load_edits(network: Network, loads: list[Load]) -> list[str]:
    """Edit commands that set the powers of the existing Load elements to loads"""
    case = network.case()
    commands = []
    for load in loads:
        names = [name for name, _, _ in dss_load_elements(case, load)]
        p, q = load.per_phase()
        edits = [
            (f"ld{load.id}_{'abc'[k]}", f"kW={p[k] / 1e3!r} kvar={q[k] / 1e3!r}")
            for k in range(3)
            if p[k] != 0.0 or q[k] != 0.0
        ]
        assert [name for name, _ in edits] == names, f"load {load.id} changes its OpenDSS elements"
        commands += [f"Edit Load.{name} {values}" for name, values in edits]
    return commands


def dss_solve() -> list[str]:
    return [
        f"Set VoltageBases={DSS_VOLTAGE_BASES}",
        "CalcVoltageBases",
        f"Set Tolerance={DSS_TOLERANCE:g}",
        f"Set MaxIterations={DSS_MAX_ITERATIONS}",
        "Solve",
    ]


def run_commands(commands: list[str]):
    for command in commands:
        dss.Text.Command(command)


# ---------------------------------------------------------------------------------------------------------------
# self-checks 1 to 4 (each element type on its own)


def element_yprim(name: str) -> np.ndarray:
    dss.Circuit.SetActiveElement(name)
    raw = np.array(dss.CktElement.YPrim())
    y = raw[0::2] + 1j * raw[1::2]
    n = math.isqrt(len(y))
    return y.reshape(n, n)


def small_circuit(kv: float, elements: list[str]):
    """a 3-phase test circuit with bus src at kv, the given elements, solved"""
    dss.Basic.ClearAll()
    run_commands(
        [
            "Clear",
            f"Set DefaultBaseFrequency={FREQ:g}",
            f"New Circuit.t phases=3 basekv={kv!r} pu=1 angle=0 bus1=src Z1=[0, 1e-6] Z0=[0, 1e-6]",
            *elements,
            "Set VoltageBases=[1]",
            "CalcVoltageBases",
            "Solve",
        ]
    )
    assert dss.Solution.Converged()


def check_linecodes(network: Network):
    """self-check 1: a 1-km line of every used LineCode. Three-phase: the Line's R, X and C matrices are the
    (2 Z1 + Z0) / 3 and (Z0 - Z1) / 3 forms (C from b / omega). One-phase: the YPrim series impedance is Z1 and the
    shunt is C1 (r0, x0, b0 ignored), which catches the semantics changing with the OpenDSS version"""
    commands = [
        f"New LineCode.lc{key} nphases={code['Phases']} r1={code['r1']!r} x1={code['x1']!r} b1={code['b1']!r}"
        f" r0={code['r0']!r} x0={code['x0']!r} b0={code['b0']!r} units=km"
        for key, code in network.linecodes.items()
    ]
    for k, (key, code) in enumerate(network.linecodes.items()):
        nodes = ".1.2.3" if code["Phases"] == 3 else ".1"
        commands.append(
            f"New Line.k{k} phases={code['Phases']} bus1=src{nodes} bus2=x{k}{nodes} linecode=lc{key} length=1 units=km"
        )
    small_circuit(22.0, commands)
    worst = {1: 0.0, 3: 0.0}
    worst_shunt = 0.0
    for k, (key, code) in enumerate(network.linecodes.items()):
        z1 = complex(code["r1"], code["x1"])
        c1 = code["b1"] * 1e-6 / OMEGA
        if code["Phases"] == 3:
            z0 = complex(code["r0"], code["x0"])
            c0 = code["b0"] * 1e-6 / OMEGA
            z_self, z_mutual = (2.0 * z1 + z0) / 3.0, (z0 - z1) / 3.0
            c_self, c_mutual = (2.0 * c1 + c0) / 3.0, (c0 - c1) / 3.0
            z_expected = np.full((3, 3), z_mutual) + np.eye(3) * (z_self - z_mutual)
            c_expected = np.full((3, 3), c_mutual) + np.eye(3) * (c_self - c_mutual)
            dss.Lines.Name(f"k{k}")
            z = (np.array(dss.Lines.RMatrix()) + 1j * np.array(dss.Lines.XMatrix())).reshape(3, 3)
            c = np.array(dss.Lines.CMatrix()).reshape(3, 3) * 1e-9
            dev = max(
                np.max(np.abs(z - z_expected)) / np.max(np.abs(z_expected)),
                np.max(np.abs(c - c_expected)) / np.max(np.abs(c_expected)),
            )
        else:
            y = element_yprim(f"Line.k{k}")
            assert y.shape == (2, 2)
            y_series = 1.0 / z1
            shunt = 0.5j * OMEGA * c1
            y_expected = np.array([[y_series + shunt, -y_series], [-y_series, y_series + shunt]])
            # the series impedance on its own, and the whole YPrim (the shunt to the round-off of y_series + shunt)
            dev = max(abs(-1.0 / y[0, 1] / z1 - 1.0), np.max(np.abs(y - y_expected)) / np.max(np.abs(y_expected)))
            worst_shunt = max(worst_shunt, abs((y[0, 0] + y[0, 1]) / shunt - 1.0))
        assert dev < 1e-12, f"LineCode {key}: OpenDSS differs from the expected line model by {dev:.1e} relative"
        worst[code["Phases"]] = max(worst[code["Phases"]], dev)
    n3 = sum(1 for code in network.linecodes.values() if code["Phases"] == 3)
    print(
        f"  check 1: {n3} three-phase LineCodes match the (2 Z1 + Z0) / 3, (Z0 - Z1) / 3 forms (max dev"
        f" {worst[3]:.1e}); {len(network.linecodes) - n3} one-phase LineCodes solve with Z1 and C1 (max dev"
        f" {worst[1]:.1e}, shunt alone {worst_shunt:.1e})"
    )


def single_phase_key(t: Transformer1P):
    return (t.center_tapped(), t.u1, t.u2, t.sn, t.uk, t.pk, t.taps())


def check_single_phase_transformers(network: Network):
    """self-check 2: each distinct ISO, SWER, phase-to-earth and split-phase parameter set on its own"""
    distinct: dict[tuple, Transformer1P] = {}
    for t in network.single_phase_transformers():
        distinct.setdefault(single_phase_key(t), t)
    for t in distinct.values():
        if t.center_tapped():
            check_center_tapped_transformer(t)
        else:
            check_transformer(t)
    print(f"  check 2: {len(distinct)} distinct single-phase transformer parameter sets agree")


def check_three_phase_transformer(t: ThreePhaseTransformer) -> tuple[float, float]:
    """self-check 3: the OpenDSS YPrim restricted to the phase conductors (the neutrals are earthed) equals the numpy
    PGM model (the decisive check of the Dy1 convention and the clock), and the mirror Reactor equals PGM's low
    susceptance -j b / 3 J (to its regularisation eps); returns the deviations"""
    elements = [f"New Transformer.t{t.id} phases=3 windings=2 buses=[src sec] {t.dss_parameters()}"]
    mirror = t.mirror_reactor("src")
    if mirror is not None:
        elements.append(mirror)
    small_circuit(t.u1 / 1e3, elements)
    y = element_yprim(f"Transformer.t{t.id}")
    assert y.shape == (8, 8)
    dss.Circuit.SetActiveElement(f"Transformer.t{t.id}")
    order = dss.CktElement.NodeOrder()
    assert list(order[0:3]) == [1, 2, 3]
    assert list(order[4:7]) == [1, 2, 3]
    phases = [0, 1, 2, 4, 5, 6]
    y_phases = y[np.ix_(phases, phases)]
    expected = t.primitive_admittance()
    dev = np.max(np.abs(y_phases - expected)) / np.max(np.abs(expected))
    assert dev < 1e-9, f"YPrim of transformer {t.id} differs from the PGM model by {dev:.1e} relative"
    dev_mirror = 0.0
    b = t.low_susceptance()
    if b is not None:
        y_mirror = element_yprim(f"Reactor.lo{t.id}src")[:3, :3]
        ones = np.ones((3, 3))
        exact = -1j * b / 3.0 * (ones + LOW_SUSCEPTANCE_EPS * np.eye(3))
        dev_exact = np.max(np.abs(y_mirror - exact)) / np.max(np.abs(exact))
        assert dev_exact < 1e-9, f"mirror Reactor of transformer {t.id}: {dev_exact:.1e} relative"
        pgm = -1j * b / 3.0 * ones
        dev_mirror = np.max(np.abs(y_mirror - pgm)) / np.max(np.abs(pgm))
        assert dev_mirror < 2.0 * LOW_SUSCEPTANCE_EPS, f"mirror Reactor of transformer {t.id}: {dev_mirror:.1e}"
    return dev, dev_mirror


def three_phase_key(t: ThreePhaseTransformer):
    return (
        t.u1,
        t.u2,
        t.sn,
        t.uk,
        t.pk,
        t.winding_from,
        t.winding_to,
        t.clock,
        t.u_rated_from,
        t.u_rated_to,
        t.tap_ratio(),
    )


def check_three_phase_transformers(network: Network):
    distinct: dict[tuple, ThreePhaseTransformer] = {}
    for t in network.three_phase_transformers():
        distinct.setdefault(three_phase_key(t), t)
    worst = 0.0
    worst_mirror = 0.0
    for t in distinct.values():
        dev, dev_mirror = check_three_phase_transformer(t)
        worst = max(worst, dev)
        worst_mirror = max(worst_mirror, dev_mirror)
    print(
        f"  check 3: {len(distinct)} distinct three-phase transformer sets (zone substation Dyn1, regulators YNyn0, LV"
        f" Dyn1): YPrim equals the PGM model (max dev {worst:.1e} relative); mirror Reactors equal PGM's low"
        f" susceptance (max dev {worst_mirror:.1e}, the eps regularisation)"
    )


def check_capacitors(network: Network):
    """self-check 4: the capacitor YPrim is j b I (the PGM shunt with b1 = b0 = b)"""
    worst = 0.0
    for c in network.capacitors:
        small_circuit(22.0, [f"New Capacitor.c{c['id']} bus1=src.1.2.3 phases=3 kvar={c['kvar']!r} kV={c['kv']!r}"])
        y = element_yprim(f"Capacitor.c{c['id']}")[:3, :3]
        expected = 1j * c["b"] * np.eye(3)
        dev = np.max(np.abs(y - expected)) / c["b"]
        assert dev < 1e-12, f"capacitor {c['id']}: YPrim differs by {dev:.1e} relative"
        worst = max(worst, dev)
    print(f"  check 4: {len(network.capacitors)} capacitors: YPrim = j b I (max dev {worst:.1e})")


# ---------------------------------------------------------------------------------------------------------------
# power flow (self-check 5) and summary


def solve(network: Network, commands: list[str], label: str) -> int:
    """run the commands (ending in Solve); abort with the D-M10-1 message if OpenDSS does not converge"""
    message = D_M10_1_MESSAGE.format(max_demand=MAX_DEMAND, network=network.spec.name)
    try:
        run_commands(commands)
    except dss.DSSException as error:
        raise SystemExit(f"{label}: OpenDSS failed ({error}): {message}") from error
    if not dss.Solution.Converged():
        raise SystemExit(f"{label}: OpenDSS did not converge in {dss.Solution.Iterations()} iterations: {message}")
    return dss.Solution.Iterations()


def check_load_powers(network: Network, loads: list[Load]) -> float:
    """self-check 5b: every Load element draws its kW + j kvar (constant power, no vminpu fallback)"""
    case = network.case()
    worst = 0.0
    for load in loads:
        p, q = load.per_phase()
        for name, _, _ in dss_load_elements(case, load):
            k = "abc".index(name[-1])
            expected = complex(p[k], q[k])
            dss.Circuit.SetActiveElement(f"Load.{name}")
            raw = np.array(dss.CktElement.Powers()).reshape(-1, 2)
            actual = complex(*raw.sum(axis=0)) * 1e3
            dev = abs(actual - expected) / abs(expected)
            assert dev < 1e-9, f"Load.{name} draws {actual} VA instead of {expected} VA"
            worst = max(worst, dev)
    return worst


def voltage_extremes(network: Network) -> dict[str, tuple[float, int, float, int]]:
    """per node class: (min pu, its node, max pu, its node) over the phases of every node"""
    names = dss.Circuit.AllNodeNames()
    magnitudes = dss.Circuit.AllBusVMag()
    result: dict[str, tuple[float, int, float, int]] = {}
    for name, magnitude in zip(names, magnitudes):
        bus, number = name.split(".")
        assert bus.startswith("n"), name
        assert number in ("1", "2", "3"), name
        node = int(bus[1:])
        node_class = network.node_class[node]
        pu = magnitude / NOMINAL_VOLTAGE[node_class]
        low, low_node, high, high_node = result.get(node_class, (math.inf, 0, -math.inf, 0))
        if pu < low:
            low, low_node = pu, node
        if pu > high:
            high, high_node = pu, node
        result[node_class] = (low, low_node, high, high_node)
    return result


def three_phase_loading(t: ThreePhaseTransformer) -> float:
    """PGM's loading of a transformer: max(sum_k |S_from,k|, sum_k |S_to,k|) / sn, with the mirror Reactor's power on
    the from side"""
    s_from, s_to = terminal_values(f"Transformer.t{t.id}", dss.CktElement.Powers, 4)
    if t.low_susceptance() is not None:
        s_low, _ = terminal_values(f"Reactor.lo{t.id}n{t.from_node}", dss.CktElement.Powers, 3)
        s_from = s_from + s_low
    return float(max(np.sum(np.abs(s_from)), np.sum(np.abs(s_to))) * 1e3 / t.sn)


def power_flow_summary(network: Network, label: str, iterations: int, loads: list[Load]):
    """self-check 5 (convergence, constant-power loads, minimum voltage per class) and the printed summary"""
    load_dev = check_load_powers(network, loads)
    extremes = voltage_extremes(network)
    total = dss.Circuit.TotalPower()
    losses = dss.Circuit.Losses()
    print(f"  {label}: converged in {iterations} iterations; loads draw their kW + j kvar (max dev {load_dev:.1e})")
    print(
        f"    feeder P = {-total[0] / 1e3:.4f} MW, Q = {-total[1] / 1e3:.4f} Mvar; losses {losses[0] / 1e3:.2f} kW;"
        f" load P = {sum(sum(load.p) for load in loads) / 1e6:.4f} MW"
    )
    for node_class in NOMINAL_VOLTAGE:
        if node_class not in extremes:
            continue
        low, low_node, high, high_node = extremes[node_class]
        print(f"    {node_class:30s} min {low:.4f} pu (node {low_node}), max {high:.4f} pu (node {high_node})")
    iso_loading = ", ".join(
        f"{network.iso_names[t.id]} {transformer_pf_reference(t)['loading']:.1%}" for t in network.isos
    )
    print(f"    ISO loading: {iso_loading}")
    regulator_loading = ", ".join(
        f"{network.regulator_names[t.id]} {three_phase_loading(t):.1%}" for t in network.regulators
    )
    print(f"    regulator loading: {regulator_loading}; zone substation {three_phase_loading(network.zone):.1%}")
    failed = [node_class for node_class, (low, *_) in extremes.items() if low < MIN_VOLTAGE_PU]
    if failed:
        message = D_M10_1_MESSAGE.format(max_demand=MAX_DEMAND, network=network.spec.name)
        raise SystemExit(f"{label}: node voltage below {MIN_VOLTAGE_PU} pu in {', '.join(failed)}: {message}")


def check_network(network: Network):
    """self-checks 1 to 5 and the dry power flow at minimum and maximum demand"""
    name = network.spec.name
    check_linecodes(network)
    check_single_phase_transformers(network)
    check_three_phase_transformers(network)
    check_capacitors(network)
    pgm_rows(network, network.loads_min)
    factors = ", ".join(
        f"{network.iso_names[iso.id]} {network.iso_factor[iso.id]:.4g}"
        f" ({iso.sn / 1e3:g} kVA for {network.iso_swer_kva[iso.id]:g} kVA)"
        for iso in network.isos
    )
    print(f"  ISO factors at maximum demand ({MAX_DEMAND:g}): {factors}")
    dss.Basic.ClearAll()
    commands = [*dss_network(network), *dss_loads(network, network.loads_min), *dss_solve()]
    iterations = solve(network, commands, f"{name} minimum demand ({MIN_DEMAND:g})")
    power_flow_summary(network, f"minimum demand ({MIN_DEMAND:g})", iterations, network.loads_min)
    commands = [*dss_load_edits(network, network.loads_max), "Solve"]
    iterations = solve(network, commands, f"{name} maximum demand ({MAX_DEMAND:g})")
    power_flow_summary(network, f"maximum demand ({MAX_DEMAND:g})", iterations, network.loads_max)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("root", type=Path, help="the tests/data directory to write to")
    parser.add_argument("networks", nargs="*", help="the networks to generate (smr8, klo14; default both)")
    parser.add_argument("--source-dir", type=Path, help="read the xlsx files from this directory")
    parser.add_argument("--check-only", action="store_true", help="run the self-checks only and write nothing")
    args = parser.parse_args()
    unknown = set(args.networks) - {spec.name for spec in NETWORKS}
    if unknown:
        parser.error(f"unknown networks {sorted(unknown)}")
    if not args.check_only:
        parser.error("writing the cases is not implemented yet (M10b); use --check-only")
    specs = [spec for spec in NETWORKS if not args.networks or spec.name in args.networks]
    print(f"opendssdirect.py {dss.__version__}; MAX_DEMAND = {MAX_DEMAND:g}, MIN_DEMAND = {MIN_DEMAND:g}")
    with tempfile.TemporaryDirectory() as download_dir:
        paths = fetch_sources(specs, args.source_dir, Path(download_dir))
        sheets = {spec.name: read_workbook(paths[spec.name]) for spec in specs}
    for spec in specs:
        network = build_network(spec, sheets[spec.name])
        check_network(network)


if __name__ == "__main__":
    main()
