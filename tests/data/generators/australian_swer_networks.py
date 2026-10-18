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
self-checks, the power flow at minimum and maximum demand and the short-circuit scenarios, prints a summary and writes
nothing.

It writes tests/data/power_flow/australian-swer/{smr8,klo14}: a batch power flow with scenario 0 at minimum and
scenario 1 at maximum demand (input.json holds the minimum demand; update_batch.json sets the asym_load powers of both
scenarios), the OpenDSS script test.dss, the reference asym_output_batch.json, params.json and a README.md, and a
.license file per file (MPL-2.0 AND BSD-3-Clause, the upstream data being BSD 3-Clause). test.dss is what was run: the
script executes it line by line and reads each reference right after the Solve of its scenario. The datasets are
written in the compact JSON format (a root "attributes" list per component and one list row per line); inputs keep
exact floats and references are rounded to 12 significant digits. The reference has node (u, u_angle), asym_line,
single_phase_transformer and transformer (p, q from and to, loading); it has no three-phase line output (the node
voltages at both ends pin every line, and it keeps the data small).

It also writes tests/data/short_circuit/australian-swer/{smr8,klo14}: a batch IEC 60909 short circuit (maximum,
c = 1.1) on the same network without loads, one fault (ID 6_000_000) re-pointed per scenario by update_batch.json. The
fault nodes are selected deterministically on the merged network (fault_scenarios; the table is in the README):
the zone-substation bus, the 22 kV node with the largest path impedance, the secondary of the largest regulator, the
SWER node with the largest path impedance from its ISO (bolted and through 20 + j5 ohm), and the LV busbars of a SWER,
the largest split-phase, the largest three-phase and (KLO14) the largest 22 kV phase-to-earth transformer. test.dss
builds the network once and, per scenario, defines the fault elements (named f6000000s{k}), solves and disables them;
the reference is read after each Solve. It has fault (i_f, i_f_angle), single_phase_transformer and transformer
(i_from, i_to and their angles), no node or line output (size). Self-checks: the Reactor of a fault with x_f carries
V / z_f, and a bolted two-phase fault solved on its own with a 10x larger bridge gives the same fault current. The same
faults are also solved on the source network, for the effect of the merge in the README.

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
- near-zero lines merged (both tools): a line with |z1| L below Z_MERGE_OHM is removed and its end nodes are merged,
  and the zone substation uk is raised to ZONE_SUBSTATION_UK; the near-zero impedances of the data otherwise make the
  admittance matrix too ill-conditioned for PGM's Newton-Raphson (see Z_MERGE_OHM and merge_connectors). The generator
  also solves the source network (every line, the sheet uk) in OpenDSS and writes the change into the README
"""

import argparse
import hashlib
import json
import math
import tempfile
import urllib.request
from collections import Counter, defaultdict, deque
from dataclasses import dataclass, field, replace
from decimal import Decimal
from pathlib import Path

import numpy as np
import opendssdirect as dss
from openpyxl import load_workbook
from reduced_phase_opendss_cases import (
    CENTER_TAPPED,
    FREQ,
    IC_ERROR_TOLERANCE,
    SC_PARAMS,
    Case,
    Fault,
    Load,
    Transformer1P,
    bus_phasors,
    bus_voltages,
    check_center_tapped_transformer,
    check_fault_impedance,
    check_transformer,
    dss_load_elements,
    fault_currents,
    fmt_matrix,
    mag_angle,
    pgm_fault,
    set_active_element,
    short_circuit_impedance,
    source_z,
    terminal_values,
    transformer_pf_reference,
    transformer_sc_reference,
)

OMEGA = 2.0 * math.pi * FREQ

# REUSE-IgnoreStart
LICENSE = (
    "SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>\n"
    "SPDX-FileCopyrightText: 2024 Team Nando - Public Repositories\n"
    "\n"
    "SPDX-License-Identifier: MPL-2.0 AND BSD-3-Clause\n"
)
UPSTREAM_COPYRIGHT = "Copyright (c) 2024, Team Nando - Public Repositories"
# REUSE-IgnoreEnd

SOURCE_REPO = "https://github.com/Team-Nando/MV-LV-Networks"
SOURCE_COMMIT = "8429f8a19f7b95acca2965cc48f2736880e5b31c"
RAW_URL = "https://raw.githubusercontent.com/Team-Nando/MV-LV-Networks/{commit}/{file}"

# demand as a fraction of the transformer kVA, at pf 0.95 lagging
POWER_FACTOR = 0.95
MIN_DEMAND = 0.10
# decision D47 (docs/advanced_documentation/multi-phase-and-swer-design.md): the 0.60 first approved is infeasible
# (OpenDSS does not converge and the voltage collapses to 0.29 pu on SMR8 and 0.53 pu on KLO14); 0.20 keeps every
# regulator within its rating
MAX_DEMAND = 0.20
# self-check 5: the lowest node voltage of every class at either demand (below it: the D47 message)
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

# Tolerances: the rule of reduced_phase_opendss_cases.py. Per case, the smallest round rtol (with atol <= 1e-6), then
# the smallest round atol, for which the worst |actual - expected| / (atol + rtol |expected|) over all methods and
# both scenarios is below 0.1 (a 10x margin), as measured by reduced_phase_compare_pgm.py. For these networks the
# power fields (p and q, from and to) get their own atol of at most 1e-3 (W, var): a branch power is computed from the
# voltages at its ends, so a short branch (down to Z_MERGE_OHM) amplifies their round-off. iterative_current runs with
# IC_ERROR_TOLERANCE (extra_params, so newton_raphson keeps its default: it stalls at about 5e-11 on these networks).
# On KLO14 it also gets its own atol of 1e-2 on the power fields (newton_raphson keeps 1e-3): the worst entry is the
# 795 W phase-c flow through the 200 kVA regulator PLENTY_VLY_RISE, where iterative_current is off by about 2e-4 W of
# round-off (newton_raphson by 1.2e-4 W), 0.11 of the 1e-3 + 1e-6 x 795 W tolerance. The cause is the fixed point of
# the iteration: it solves the rounded matrix, so the LU round-off is not removed by iterating (err_tol 1e-9 to 1e-14
# all give the same 0.11), whereas newton_raphson zeroes a residual that uses the true matrix. A defect-correction
# (residual) form of iterative_current was tried and rejected: it only turned it into a residual-noise-limited method
# like newton_raphson (stall at 3e-11 to 5e-11, ratios moving with the tolerance) without a systematic gain. With the
# IC-only atol the worst ratio of iterative_current is 0.018. SMR8 needs no IC-only atol; its IC_ERROR_TOLERANCE run
# agrees as well as newton_raphson (worst ratios 0.028 and 0.050), which tightens its rtol from 1e-6 (the rtol that
# iterative_current needed at the default error tolerance 1e-8, ratio 0.066) to 1e-7. KLO14 stays at 1e-6, which
# newton_raphson needs (rtol 1e-7 gives 0.137).
POWER_FIELDS = "[pq]_(from|to)"
KLO14_IC_POWER_ATOL = 1e-2
IC_EXTRA_PARAMS = {"error_tolerance": IC_ERROR_TOLERANCE}
PF_PARAMS = {
    "smr8": {
        "calculation_method": ["newton_raphson", "iterative_current"],
        "rtol": 1e-7,
        "atol": {"default": 1e-8, POWER_FIELDS: 1e-3},
        "experimental_features": "enabled",
        "extra_params": {"iterative_current": dict(IC_EXTRA_PARAMS)},
    },
    "klo14": {
        "calculation_method": ["newton_raphson", "iterative_current"],
        "rtol": 1e-6,
        "atol": {"default": 1e-8, POWER_FIELDS: 1e-3},
        "experimental_features": "enabled",
        "extra_params": {
            "iterative_current": dict(IC_EXTRA_PARAMS, atol={"default": 1e-8, POWER_FIELDS: KLO14_IC_POWER_ATOL})
        },
    },
}
# the measured agreement, for the README
PF_AGREEMENT = {
    "smr8": (
        "Measured: newton_raphson converges in 3 and 4 iterations and iterative_current (error_tolerance 1e-10) in 10"
        " and 14 (minimum, maximum demand); the node voltages agree to 1.1e-10 (newton_raphson) and 9.6e-11"
        " (iterative_current) relative and the powers to 1.0e-3 W absolute, worst ratios of 0.050 (newton_raphson)"
        " and 0.044 (iterative_current). At the default error tolerance 1e-8 iterative_current needed rtol 1e-6"
        " (ratio 0.066)."
    ),
    "klo14": (
        "Measured: newton_raphson converges in 3 and 4 iterations (minimum, maximum demand); the node voltages agree"
        " to 9.8e-11 relative and the powers to 1.2e-3 W absolute, a worst ratio of 0.069. iterative_current (error"
        " tolerance 1e-10, 9 and 12 iterations) has its own atol of 1e-2 on the power fields, set in `extra_params`"
        " so that newton_raphson keeps 1e-3 (worst ratios 0.018 and 0.016 for iterative_current): its worst entry is"
        " the 795 W phase-c flow through the 200 kVA regulator PLENTY_VLY_RISE, where it is off by about 2e-4 W of"
        " round-off (newton_raphson 1.2e-4 W), a ratio of 0.11 against 1e-3. The round-off comes from the fixed point"
        " of the iteration, which solves the rounded matrix, so a tighter error tolerance does not remove it. A"
        " residual (defect-correction) form of iterative_current was tried and rejected: it did not agree better with"
        " newton_raphson."
    ),
}

# Short circuit: IEC 60909 maximum (c = 1.1, as PGM takes it from the 66 kV source node), no loads. The fault
# impedances of the scenarios: a two-phase fault on 22 kV gets r_f = 0.5 ohm (a bolted 22 kV bridge has a round-off of
# about 1.6e-2 in OpenDSS, see reduced_phase_opendss_cases.py), and the high-impedance SWER fault r_f + j x_f.
SC_VOLTAGE_FACTOR = 1.1
TWO_PHASE_R_F = 0.5
SWER_FAULT_Z = complex(20.0, 5.0)
# Tolerances, by the rule above: rtol 1e-9 misses it (the fault current angles agree to about 3e-10 relative). The
# branch currents (i_from, i_to and their angles) get their own atol of 1e-6 A: a small current through a
# low-impedance branch (the charging currents through the ISOs and regulators, the 66 kV current of the zone
# substation for a remote fault) is computed from a small voltage difference, which amplifies the voltage round-off of
# about 1e-10 to a few 1e-8 A.
BRANCH_CURRENT_FIELDS = "i_(from|to)(_angle)?"
SC_CASE_PARAMS = {
    name: dict(SC_PARAMS, rtol=1e-8, atol={"default": 1e-8, BRANCH_CURRENT_FIELDS: 1e-6}) for name in ("smr8", "klo14")
}
# the measured agreement, for the README
SC_AGREEMENT = {
    "smr8": (
        "Measured: the fault current magnitudes agree to 6.7e-11 relative and the phasors to 3.0e-10, the branch"
        " current magnitudes above 1 A to 5.4e-9 relative and the currents up to 10 A to 3.6e-8 A absolute (as"
        " phasors), a worst ratio of 0.033 (the 66 kV current of the zone substation in scenario 4, 9.7 A)."
    ),
    "klo14": (
        "Measured: the fault current magnitudes agree to 9.2e-11 relative and the phasors to 2.9e-10, the branch"
        " current magnitudes above 1 A to 6.0e-9 relative and the currents up to 10 A to 3.2e-8 A absolute (as"
        " phasors), a worst ratio of 0.030 (the 66 kV current of the zone substation in scenario 3, 3.8 A)."
    ),
}

# Merging of near-zero lines. The source data has near-zero impedances: connector lines down to 1.4e-7 ohm
# (three-phase) and 1.4e-6 ohm (one-phase), and a zone substation with uk = 1.41e-6. They make the admittance matrix so
# ill-conditioned that PGM's Newton-Raphson stalls at a max deviation of about 3e-7 (error tolerance 1e-8), also on the
# three-phase part alone (the upstream code path), and PGM agrees with OpenDSS to only about 2e-7. A line whose series
# impedance |z1| L is below Z_MERGE_OHM is removed and its two end nodes are merged into one (in both tools); see
# merge_connectors. The zone substation gets ZONE_SUBSTATION_UK (r / x kept) if its sheet uk is smaller.
# Measured (OpenDSS, against the source network; PGM newton_raphson with the default error tolerance 1e-8):
# - Z_MERGE_OHM 1e-4 changes the node voltages by at most 1.0e-5 relative, but the tolerance rule fails on both
#   networks (the smallest line left is 1.5e-4 ohm), also with ZONE_SUBSTATION_UK 1e-4; 1e-3 changes them by 5.5e-5 and
#   the losses by up to 0.1 %, and passes the rule with ZONE_SUBSTATION_UK 1e-4; 1e-2 changes them by 1.6e-3 and the
#   losses by up to 3.2 %.
# - with uk 1e-5 (or 3e-5) the zone substation (5.2e-5 ohm) is the smallest impedance left after the merge, and KLO14
#   misses the rule (worst ratio 0.145 at Z_MERGE_OHM 1e-3, the powers through a 200 kVA regulator); uk 1e-4
#   (5.2e-4 ohm) passes it, and on its own changes the node voltages by at most 1.3e-5 and the losses by 0.013 %.
# Newton-Raphson converges in 3 and 4 iterations with every one of these values.
Z_MERGE_OHM = 1e-3
ZONE_SUBSTATION_UK = 1e-4

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

D47_MESSAGE = (
    "maximum demand MAX_DEMAND = {max_demand} of the transformer kVA is infeasible on {network} (decision D47 in"
    " docs/advanced_documentation/multi-phase-and-swer-design.md: 0.60 does not converge in OpenDSS and collapses the"
    " voltage; 0.20 is the chosen maximum)"
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


@dataclass(frozen=True)
class MergeInfo:
    """what merge_connectors changed"""

    representative: dict[int, int]  # every merged node -> the node it is merged into
    merged_lines: dict[int, int]  # the number of removed (merged) lines by phases
    dropped_loops: dict[int, int]  # the number of lines that became a self-loop by phases
    parallel_lines: int  # lines that run in parallel to another line after the merge
    largest_merged: float  # ohm, the largest |z1| L of a merged line
    smallest_kept: float  # ohm, the smallest |z1| L of a kept line
    nodes_before: int
    lines_before: dict[int, int]  # by phases
    zone_sheet_uk: float


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
    merge: MergeInfo | None = None  # None for the source network

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


def line_impedance(line: LineSection, linecodes: dict[str, dict]) -> float:
    """the series impedance |z1| L of a line (ohm)"""
    code = linecodes[line.code]
    z_per_km = math.hypot(code["r1"], code["x1"])
    assert z_per_km > 0.0, f"LineCode {line.code} has no series impedance"
    return z_per_km * line.length


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
    """the source network: the section 2 mapping of the plan, with the topology checks (self-check 7), before the merge
    of near-zero lines"""
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
                    # decision D52: both data ratings taken phase-to-earth (the data's turns ratio)
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

    # demand (decisions D47 and D48)
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


def network_edges(network: Network) -> list[tuple[int, int, int]]:
    """every branch as (id, from node, to node): lines and transformers"""
    return [
        *((line.id, line.from_node, line.to_node) for line in network.lines),
        *((t.id, t.from_node, t.to_node) for t in network.three_phase_transformers()),
        *((t.id, t.from_node, t.to_node) for t in network.single_phase_transformers()),
    ]


def source_distances(network: Network) -> dict[int, int]:
    """the number of branches between the source bus and every node (breadth first); every node is reached"""
    adjacency = defaultdict(list)
    for _, a, b in network_edges(network):
        adjacency[a].append(b)
        adjacency[b].append(a)
    distance = {SOURCE_NODE: 0}
    queue = deque([SOURCE_NODE])
    while queue:
        node = queue.popleft()
        for other in adjacency[node]:
            if other not in distance:
                distance[other] = distance[node] + 1
                queue.append(other)
    assert distance.keys() == network.node_u_rated.keys(), "a node is not reachable from the source"
    return distance


def rebus(bus: str, old: int, new: int) -> str:
    """an OpenDSS bus n{old}.x.y as n{new}.x.y"""
    prefix = f"n{old}."
    assert bus.startswith(prefix), bus
    return f"n{new}.{bus[len(prefix) :]}"


def merge_connectors(source: Network) -> Network:
    """the network as solved by both tools: every line with |z1| L below Z_MERGE_OHM is removed and its two end nodes
    are merged into one (union-find). A merged group of nodes is represented by its node closest to the source (the
    fewest branches; on a tie the lowest id), and every element on a merged node moves to the representative: lines,
    transformer terminals (zone substation, regulators, ISOs, distribution transformers) and capacitors. A kept line
    whose ends merge is dropped (a self-loop: it runs in parallel with merged lines only), and lines that end up
    between the same pair of nodes stay parallel. A three-phase line merges 22 kV nodes and a one-phase line SWER
    nodes; a 22 kV node never merges with a SWER node. The zone substation uk is raised to ZONE_SUBSTATION_UK if its
    sheet uk is smaller (uk and pk scaled together, r / x kept)."""
    distance = source_distances(source)
    parent = {node: node for node in source.node_u_rated}

    def find(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    expected_class = {3: CLASS_MV, 1: CLASS_SWER}
    merged_ids = set()
    for line in source.lines:
        if line_impedance(line, source.linecodes) >= Z_MERGE_OHM:
            continue
        ends = {source.node_class[line.from_node], source.node_class[line.to_node]}
        assert ends == {expected_class[line.phases]}, f"line {line.id} ({line.phases} phases) between {ends}"
        merged_ids.add(line.id)
        a, b = find(line.from_node), find(line.to_node)
        if a != b:
            keep, drop = sorted((a, b), key=lambda node: (distance[node], node))
            parent[drop] = keep
    representative = {node: find(node) for node in parent if find(node) != node}
    for node, kept in representative.items():
        assert source.node_class[node] == source.node_class[kept], (
            f"node {node} ({source.node_class[node]}) merges with node {kept} ({source.node_class[kept]})"
        )
        assert source.node_class[node] in (CLASS_MV, CLASS_SWER), node

    def rep(node: int) -> int:
        return representative.get(node, node)

    lines = []
    dropped_loops: Counter[int] = Counter()
    for line in source.lines:
        if line.id in merged_ids:
            continue
        a, b = rep(line.from_node), rep(line.to_node)
        if a == b:
            dropped_loops[line.phases] += 1
        else:
            lines.append(replace(line, from_node=a, to_node=b))
    pairs = Counter(frozenset((line.from_node, line.to_node)) for line in lines)
    parallel_lines = sum(n for n in pairs.values() if n > 1)

    zone_uk, zone_pk = source.zone.uk, source.zone.pk
    if zone_uk < ZONE_SUBSTATION_UK:
        zone_pk *= ZONE_SUBSTATION_UK / zone_uk
        zone_uk = ZONE_SUBSTATION_UK
    zone = replace(source.zone, to_node=rep(source.zone.to_node), uk=zone_uk, pk=zone_pk)
    regulators = [replace(t, from_node=rep(t.from_node), to_node=rep(t.to_node)) for t in source.regulators]
    isos = [replace(t, from_node=rep(t.from_node), to_node=rep(t.to_node)) for t in source.isos]
    iso_buses = {t.id: rebus(source.iso_buses[t.id], t.from_node, rep(t.from_node)) for t in source.isos}
    distribution = [replace(t, from_node=rep(t.from_node)) for t in source.distribution]
    distribution_buses = {
        t.id: [rebus(buses[0], t.from_node, rep(t.from_node)), *buses[1:]]
        for t in source.distribution
        if (buses := source.distribution_buses.get(t.id)) is not None
    }
    assert all(t.from_node != t.to_node for t in [zone, *regulators, *isos, *distribution])
    assert all(load.node not in representative for load in source.loads_min)
    used_codes = {line.code for line in lines}
    network = replace(
        source,
        linecodes={key: code for key, code in source.linecodes.items() if key in used_codes},
        lines=lines,
        node_class={node: c for node, c in source.node_class.items() if node not in representative},
        node_u_rated={node: u for node, u in source.node_u_rated.items() if node not in representative},
        zone=zone,
        regulators=regulators,
        isos=isos,
        iso_buses=iso_buses,
        distribution=distribution,
        distribution_buses=distribution_buses,
        capacitors=[{**c, "node": rep(c["node"])} for c in source.capacitors],
        merge=MergeInfo(
            representative=representative,
            merged_lines=dict(Counter(line.phases for line in source.lines if line.id in merged_ids)),
            dropped_loops=dict(dropped_loops),
            parallel_lines=parallel_lines,
            largest_merged=max(
                line_impedance(line, source.linecodes) for line in source.lines if line.id in merged_ids
            ),
            smallest_kept=min(line_impedance(line, source.linecodes) for line in lines),
            nodes_before=len(source.node_u_rated),
            lines_before=dict(Counter(line.phases for line in source.lines)),
            zone_sheet_uk=source.zone.uk,
        ),
    )
    # the merged network is connected, and every node is used
    source_distances(network)
    check_unique_pgm_ids(network)
    print(f"  merge: {merge_text(network)}")
    return network


def merge_text(network: Network) -> str:
    """what merge_connectors changed, for the summary and the README"""
    merge = network.merge
    assert merge is not None
    after = Counter(line.phases for line in network.lines)
    loops = sum(merge.dropped_loops.values())
    zone_z = network.zone.uk * network.zone.u2**2 / network.zone.sn
    if network.zone.uk != merge.zone_sheet_uk:
        zone = (
            f"The zone substation uk is raised from {merge.zone_sheet_uk:.3g} to {network.zone.uk:g} ({zone_z:.2g} ohm"
            f" on the {network.zone.u2 / 1e3:g} kV side)."
        )
    else:
        zone = f"The zone substation keeps its sheet uk of {network.zone.uk:.3g}."
    return (
        f"{merge.merged_lines.get(3, 0)} three-phase and {merge.merged_lines.get(1, 0)} one-phase lines are merged"
        f" (|z1| L up to {merge.largest_merged:.3g} ohm), and {len(merge.representative)} nodes with them, which"
        f" leaves {len(network.node_u_rated)} of {merge.nodes_before} nodes and {after[3]} of"
        f" {merge.lines_before.get(3, 0)} three-phase and {after[1]} of {merge.lines_before.get(1, 0)} one-phase lines;"
        f" {loops} kept line{'s' if loops != 1 else ''} became a self-loop and {'were' if loops != 1 else 'was'}"
        f" dropped, and {merge.parallel_lines} lines run in parallel after the merge. The smallest |z1| L left is"
        f" {merge.smallest_kept:.3g} ohm. {zone}"
    )


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


def pgm_line(line: LineSection, code: dict) -> dict[str, int | float]:
    """the PGM row of a line (b of the LineCode in uS/km): a three-phase line is a line with the sequence values of its
    LineCode, a one-phase line an asym_line on phase a with Z1 and C1 only (as OpenDSS solves a 1-phase LineCode)"""
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
    else:
        row.update(
            {
                "phases": 1,
                "r_aa": code["r1"] * length,
                "x_aa": code["x1"] * length,
                "c_aa": code["b1"] * 1e-6 / OMEGA * length,
                "i_n": i_n,
            }
        )
    return row


def pgm_rows(network: Network, loads: list[Load]) -> dict[str, list[dict]]:
    """the PGM input of the network with the given loads, per component, sorted by id"""
    case = network.case()
    rows: dict[int, list[dict]] = {3: [], 1: []}
    for line in network.lines:
        rows[line.phases].append(pgm_line(line, network.linecodes[line.code]))
    lines, asym_lines = rows[3], rows[1]
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


def dss_linecode(key: str, code: dict) -> str:
    """the OpenDSS LineCode lc{key}, with the sheet values verbatim (b in uS/km)"""
    return (
        f"New LineCode.lc{key} nphases={code['Phases']} r1={code['r1']!r} x1={code['x1']!r} b1={code['b1']!r}"
        f" r0={code['r0']!r} x0={code['x0']!r} b0={code['b0']!r} units=km"
    )


def dss_line(line: LineSection, from_bus: str, to_bus: str) -> str:
    """the OpenDSS Line l{id} on its LineCode, on phases 1, 2, 3 or (one-phase) 1 of the buses"""
    nodes = ".1.2.3" if line.phases == 3 else ".1"
    return (
        f"New Line.l{line.id} phases={line.phases} bus1={from_bus}{nodes} bus2={to_bus}{nodes}"
        f" linecode=lc{line.code} length={line.length!r} units=km"
    )


def dss_capacitor(capacitor: dict, bus: str) -> str:
    """the OpenDSS Capacitor c{id}: three-phase, grounded wye (the default)"""
    return (
        f"New Capacitor.c{capacitor['id']} bus1={bus}.1.2.3 phases=3 kvar={capacitor['kvar']!r} kV={capacitor['kv']!r}"
    )


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
    commands += [dss_linecode(key, code) for key, code in network.linecodes.items()]
    commands += [dss_line(line, f"n{line.from_node}", f"n{line.to_node}") for line in network.lines]
    commands += [dss_capacitor(c, f"n{c['node']}") for c in sorted(network.capacitors, key=lambda c: c["id"])]
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
    set_active_element(name)
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
    """self-check 1: a 1-km line of every used LineCode, written as in the network, against the PGM row of the same
    line (pgm_line). Three-phase: the Line's R, X and C matrices are the (2 Z1 + Z0) / 3 and (Z0 - Z1) / 3 forms of
    the row's r1, x1, c1, r0, x0, c0. One-phase: the YPrim series impedance is r_aa + j x_aa and the shunt c_aa (r0,
    x0, b0 ignored), which catches the semantics changing with the OpenDSS version"""
    test_lines = [
        LineSection(id=k, from_node=0, to_node=0, phases=code["Phases"], length=1.0, code=key)
        for k, (key, code) in enumerate(network.linecodes.items())
    ]
    small_circuit(
        22.0,
        [
            *(dss_linecode(key, code) for key, code in network.linecodes.items()),
            *(dss_line(line, "src", f"x{line.id}") for line in test_lines),
        ],
    )
    worst = {1: 0.0, 3: 0.0}
    worst_shunt = 0.0
    for line in test_lines:
        key = line.code
        row = pgm_line(line, network.linecodes[key])
        if line.phases == 3:
            z1, z0 = complex(row["r1"], row["x1"]), complex(row["r0"], row["x0"])
            c1, c0 = row["c1"], row["c0"]
            z_self, z_mutual = (2.0 * z1 + z0) / 3.0, (z0 - z1) / 3.0
            c_self, c_mutual = (2.0 * c1 + c0) / 3.0, (c0 - c1) / 3.0
            z_expected = np.full((3, 3), z_mutual) + np.eye(3) * (z_self - z_mutual)
            c_expected = np.full((3, 3), c_mutual) + np.eye(3) * (c_self - c_mutual)
            dss.Lines.Name(f"l{line.id}")
            z = (np.array(dss.Lines.RMatrix()) + 1j * np.array(dss.Lines.XMatrix())).reshape(3, 3)
            c = np.array(dss.Lines.CMatrix()).reshape(3, 3) * 1e-9
            dev = max(
                np.max(np.abs(z - z_expected)) / np.max(np.abs(z_expected)),
                np.max(np.abs(c - c_expected)) / np.max(np.abs(c_expected)),
            )
        else:
            z1, c1 = complex(row["r_aa"], row["x_aa"]), row["c_aa"]
            y = element_yprim(f"Line.l{line.id}")
            assert y.shape == (2, 2)
            y_series = 1.0 / z1
            shunt = 0.5j * OMEGA * c1
            y_expected = np.array([[y_series + shunt, -y_series], [-y_series, y_series + shunt]])
            # the series impedance on its own, and the whole YPrim (the shunt to the round-off of y_series + shunt)
            dev = max(abs(-1.0 / y[0, 1] / z1 - 1.0), np.max(np.abs(y - y_expected)) / np.max(np.abs(y_expected)))
            worst_shunt = max(worst_shunt, abs((y[0, 0] + y[0, 1]) / shunt - 1.0))
        assert dev < 1e-12, f"LineCode {key}: OpenDSS differs from the PGM line by {dev:.1e} relative"
        worst[line.phases] = max(worst[line.phases], dev)
    n3 = sum(1 for code in network.linecodes.values() if code["Phases"] == 3)
    print(
        f"  check 1: {n3} three-phase LineCodes match the (2 Z1 + Z0) / 3, (Z0 - Z1) / 3 forms of the PGM line (max dev"
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
    set_active_element(f"Transformer.t{t.id}")
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
        small_circuit(22.0, [dss_capacitor(c, "src")])
        y = element_yprim(f"Capacitor.c{c['id']}")[:3, :3]
        expected = 1j * c["b"] * np.eye(3)
        dev = np.max(np.abs(y - expected)) / c["b"]
        assert dev < 1e-12, f"capacitor {c['id']}: YPrim differs by {dev:.1e} relative"
        worst = max(worst, dev)
    print(f"  check 4: {len(network.capacitors)} capacitors: YPrim = j b I (max dev {worst:.1e})")


# ---------------------------------------------------------------------------------------------------------------
# power flow (self-check 5) and summary


def solve(network: Network, commands: list[str], label: str) -> int:
    """run the commands (ending in Solve); abort with the D47 message if OpenDSS does not converge"""
    message = D47_MESSAGE.format(max_demand=MAX_DEMAND, network=network.spec.name)
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
            set_active_element(f"Load.{name}")
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


def three_phase_pf_reference(t: ThreePhaseTransformer) -> dict:
    """per-phase powers of both sides of a PGM transformer (the neutral conductors dropped), with the mirror Reactor's
    power on the from side, and PGM's loading max(sum_k |S_from,k|, sum_k |S_to,k|) / sn"""
    s_from, s_to = terminal_values(f"Transformer.t{t.id}", dss.CktElement.Powers, 4)
    if t.low_susceptance() is not None:
        s_low, _ = terminal_values(f"Reactor.lo{t.id}n{t.from_node}", dss.CktElement.Powers, 3)
        s_from = s_from + s_low
    s_from = s_from * 1e3
    s_to = s_to * 1e3
    return {
        "id": t.id,
        "p_from": [float(v.real) for v in s_from],
        "q_from": [float(v.imag) for v in s_from],
        "p_to": [float(v.real) for v in s_to],
        "q_to": [float(v.imag) for v in s_to],
        "loading": float(max(np.sum(np.abs(s_from)), np.sum(np.abs(s_to))) / t.sn),
    }


def three_phase_loading(t: ThreePhaseTransformer) -> float:
    return three_phase_pf_reference(t)["loading"]


@dataclass
class PowerFlowSummary:
    label: str
    iterations: int
    feeder_p: float  # W
    feeder_q: float  # var
    losses: float  # W
    load_p: float  # W
    extremes: dict[str, tuple[float, int, float, int]]
    iso_loading: dict[int, float]
    regulator_loading: dict[int, float]
    zone_loading: float


def power_flow_summary(network: Network, label: str, iterations: int, loads: list[Load]) -> PowerFlowSummary:
    """self-check 5 (convergence, constant-power loads, minimum voltage per class) and the printed summary"""
    load_dev = check_load_powers(network, loads)
    total = dss.Circuit.TotalPower()
    losses = dss.Circuit.Losses()
    summary = PowerFlowSummary(
        label=label,
        iterations=iterations,
        feeder_p=-total[0] * 1e3,
        feeder_q=-total[1] * 1e3,
        losses=losses[0],
        load_p=sum(sum(load.p) for load in loads),
        extremes=voltage_extremes(network),
        iso_loading={t.id: transformer_pf_reference(t)["loading"] for t in network.isos},
        regulator_loading={t.id: three_phase_loading(t) for t in network.regulators},
        zone_loading=three_phase_loading(network.zone),
    )
    print(f"  {label}: converged in {iterations} iterations; loads draw their kW + j kvar (max dev {load_dev:.1e})")
    print(
        f"    feeder P = {summary.feeder_p / 1e6:.4f} MW, Q = {summary.feeder_q / 1e6:.4f} Mvar; losses"
        f" {summary.losses / 1e3:.2f} kW; load P = {summary.load_p / 1e6:.4f} MW"
    )
    for node_class in NOMINAL_VOLTAGE:
        if node_class not in summary.extremes:
            continue
        low, low_node, high, high_node = summary.extremes[node_class]
        print(f"    {node_class:30s} min {low:.4f} pu (node {low_node}), max {high:.4f} pu (node {high_node})")
    iso_loading = ", ".join(f"{network.iso_names[i]} {v:.1%}" for i, v in summary.iso_loading.items())
    print(f"    ISO loading: {iso_loading}")
    regulator_loading = ", ".join(f"{network.regulator_names[i]} {v:.1%}" for i, v in summary.regulator_loading.items())
    print(f"    regulator loading: {regulator_loading}; zone substation {summary.zone_loading:.1%}")
    failed = [node_class for node_class, (low, *_) in summary.extremes.items() if low < MIN_VOLTAGE_PU]
    if failed:
        message = D47_MESSAGE.format(max_demand=MAX_DEMAND, network=network.spec.name)
        raise SystemExit(f"{label}: node voltage below {MIN_VOLTAGE_PU} pu in {', '.join(failed)}: {message}")
    return summary


# ---------------------------------------------------------------------------------------------------------------
# reference output


def rounded(value):
    """a reference value rounded to 12 significant digits (an error of at most 5e-13 relative)"""
    if isinstance(value, float):
        return float(f"{value:.12g}")
    if isinstance(value, list):
        return [rounded(v) for v in value]
    if isinstance(value, dict):
        return {k: rounded(v) for k, v in value.items()}
    return value


def pf_reference(network: Network) -> dict[str, list[dict]]:
    """the reference of the solved power flow: node, asym_line, single_phase_transformer and transformer (no
    three-phase line), each sorted by id and rounded to 12 significant digits"""
    nodes = []
    for node in sorted(network.node_u_rated):
        u, angle = bus_voltages(node)
        nodes.append({"id": node, "u": u, "u_angle": angle})
    asym_lines = []
    for line in network.lines:
        if line.phases != 1:
            continue
        s_from, s_to = terminal_values(f"Line.l{line.id}", dss.CktElement.Powers, 1)
        s_from = s_from * 1e3
        s_to = s_to * 1e3
        asym_lines.append(
            {
                "id": line.id,
                "p_from": [float(v.real) for v in s_from],
                "q_from": [float(v.imag) for v in s_from],
                "p_to": [float(v.real) for v in s_to],
                "q_to": [float(v.imag) for v in s_to],
            }
        )
    return rounded(
        {
            "node": nodes,
            "asym_line": asym_lines,
            "single_phase_transformer": [
                transformer_pf_reference(t) for t in sorted(network.single_phase_transformers(), key=lambda t: t.id)
            ],
            "transformer": [
                three_phase_pf_reference(t) for t in sorted(network.three_phase_transformers(), key=lambda t: t.id)
            ],
        }
    )


# ---------------------------------------------------------------------------------------------------------------
# self-checks and the power flow


def check_elements(network: Network):
    """self-checks 1 to 4 and the ISO factors"""
    check_linecodes(network)
    check_single_phase_transformers(network)
    check_three_phase_transformers(network)
    check_capacitors(network)
    factors = ", ".join(
        f"{network.iso_names[iso.id]} {network.iso_factor[iso.id]:.4g}"
        f" ({iso.sn / 1e3:g} kVA for {network.iso_swer_kva[iso.id]:g} kVA)"
        for iso in network.isos
    )
    print(f"  ISO factors at maximum demand ({MAX_DEMAND:g}): {factors}")


@dataclass
class PowerFlowRun:
    script: list[str]  # test.dss: exactly the commands that were run
    references: list[dict]  # per scenario
    summaries: list[PowerFlowSummary]  # per scenario


@dataclass
class SourceDeviation:
    """the OpenDSS solution of the merged network against the source network (every line of the data, the sheet uk
    of the zone substation), per scenario"""

    voltage: list[dict[str, float]]  # per node class: the max |V - V_source| / |V_source| over the kept nodes
    voltage_merged_nodes: list[float]  # the same over the merged nodes (V at their representative)
    losses: list[tuple[float, float]]  # W: (source, merged)
    feeder_p: list[tuple[float, float]]  # W: (source, merged)

    def max_voltage(self, k: int) -> float:
        return max(self.voltage[k].values())


def dss_states(network: Network) -> list[tuple[dict[str, complex], float, float]]:
    """per scenario, the OpenDSS node voltages by node name, the losses and the feeder P (W)"""
    dss.Basic.ClearAll()
    states = []
    for commands in (
        [*dss_network(network), *dss_loads(network, network.loads_min), *dss_solve()],
        [*dss_load_edits(network, network.loads_max), "Solve"],
    ):
        run_commands(commands)
        assert dss.Solution.Converged()
        raw = np.array(dss.Circuit.AllBusVolts())
        voltages = dict(zip(dss.Circuit.AllNodeNames(), raw[0::2] + 1j * raw[1::2]))
        states.append((voltages, dss.Circuit.Losses()[0], -dss.Circuit.TotalPower()[0] * 1e3))
    return states


def source_deviation(network: Network, source: Network) -> SourceDeviation:
    """how much the merge (and the zone substation uk) changes the OpenDSS solution against the source network: the
    node voltage phasors (relative, per class over the kept nodes, and over the merged nodes against their
    representative), the losses and the feeder P"""
    merge = network.merge
    assert merge is not None
    deviation = SourceDeviation(voltage=[], voltage_merged_nodes=[], losses=[], feeder_p=[])
    for (v_source, loss_source, p_source), (v_merged, loss_merged, p_merged) in zip(
        dss_states(source), dss_states(network)
    ):
        per_class: dict[str, float] = {}
        merged_nodes = 0.0
        for name, value in v_source.items():
            bus, number = name.split(".")
            node = int(bus[1:])
            rep = merge.representative.get(node, node)
            dev = float(abs(v_merged[f"n{rep}.{number}"] - value) / abs(value))
            if rep == node:
                node_class = source.node_class[node]
                per_class[node_class] = max(per_class.get(node_class, 0.0), dev)
            else:
                merged_nodes = max(merged_nodes, dev)
        deviation.voltage.append(per_class)
        deviation.voltage_merged_nodes.append(merged_nodes)
        deviation.losses.append((loss_source, loss_merged))
        deviation.feeder_p.append((p_source, p_merged))
    for k, kind in enumerate(("minimum", "maximum")):
        (loss_source, loss_merged), (p_source, p_merged) = deviation.losses[k], deviation.feeder_p[k]
        print(
            f"  against the source network at {kind} demand: node voltages {deviation.max_voltage(k):.2e}"
            f" (merged nodes {deviation.voltage_merged_nodes[k]:.2e}); losses {loss_source / 1e3:.3f} ->"
            f" {loss_merged / 1e3:.3f} kW ({loss_merged / loss_source - 1.0:+.2e}); feeder P"
            f" {p_merged / p_source - 1.0:+.2e}"
        )
    return deviation


def deviation_table(deviation: SourceDeviation) -> list[str]:
    lines = [
        "| change against the source network | scenario 0: minimum demand | scenario 1: maximum demand |",
        "|---|---|---|",
    ]
    for node_class in NOMINAL_VOLTAGE:
        if node_class in deviation.voltage[0]:
            values = " | ".join(f"{v[node_class]:.1e}" for v in deviation.voltage)
            lines.append(f"| node voltage, {node_class} | {values} |")
    values = " | ".join(f"{v:.1e}" for v in deviation.voltage_merged_nodes)
    lines.append(f"| node voltage, merged nodes | {values} |")
    values = " | ".join(
        f"{source / 1e3:.3f} to {merged / 1e3:.3f} kW ({merged / source - 1.0:+.1e})"
        for source, merged in deviation.losses
    )
    lines.append(f"| losses | {values} |")
    values = " | ".join(f"{merged / source - 1.0:+.1e}" for source, merged in deviation.feeder_p)
    lines.append(f"| feeder P (relative) | {values} |")
    return lines


def run_power_flow(network: Network) -> PowerFlowRun:
    """the batch power flow, as test.dss runs it: the network and the minimum-demand Loads, Solve (scenario 0), then
    Edit of every Load to its maximum demand, Solve (scenario 1); each reference is read right after its Solve, and
    self-check 5 runs on each scenario"""
    name = network.spec.name
    dss.Basic.ClearAll()
    run = PowerFlowRun(
        script=[
            (
                f"! power flow on the Australian SWER network {name}; written by tests/data/generators/"
                "australian_swer_networks.py, which runs these commands line by line"
            ),
            "! scenario 0 is minimum demand, scenario 1 maximum demand (see README.md)",
        ],
        references=[],
        summaries=[],
    )
    scenarios = ((MIN_DEMAND, "minimum", network.loads_min), (MAX_DEMAND, "maximum", network.loads_max))
    for k, (fraction, kind, loads) in enumerate(scenarios):
        label = f"{kind} demand ({fraction:g})"
        if k == 0:
            commands = [*dss_network(network), f"! scenario {k}: {label}", *dss_loads(network, loads), *dss_solve()]
        else:
            commands = [f"! scenario {k}: {label}", *dss_load_edits(network, loads), "Solve"]
        if k == 0:
            run_commands(run.script)
        iterations = solve(network, commands, f"{name} {label}")
        run.script += commands
        run.summaries.append(power_flow_summary(network, label, iterations, loads))
        run.references.append(pf_reference(network))
    return run


# ---------------------------------------------------------------------------------------------------------------
# short circuit


@dataclass
class ScenarioFault(Fault):
    """a Fault whose OpenDSS elements are named per scenario (f{id}s{scenario}), so that every scenario of test.dss
    defines new elements and disables them after its Solve; the PGM fault ID stays FAULT_ID"""

    scenario: int = 0

    def dss_elements(self, bridge_scale=1.0):
        prefix = f".f{self.id}"
        return [
            (element.replace(prefix, f"{prefix}s{self.scenario}", 1), definition)
            for element, definition in super().dss_elements(bridge_scale)
        ]


@dataclass(frozen=True)
class FaultScenario:
    fault: ScenarioFault
    location: str  # where, for the README
    reason: str  # why this node, for the README


def path_impedances(network: Network) -> dict[int, float]:
    """the series impedance sum |z1| L of the lines on the path from the source bus to every node, along the breadth
    first tree of the network (the branches of a node in the order of their ids; a transformer counts 0)"""
    z_line = {line.id: line_impedance(line, network.linecodes) for line in network.lines}
    adjacency = defaultdict(list)
    for branch_id, a, b in sorted(network_edges(network)):
        adjacency[a].append((branch_id, b))
        adjacency[b].append((branch_id, a))
    path = {SOURCE_NODE: 0.0}
    queue = deque([SOURCE_NODE])
    while queue:
        node = queue.popleft()
        for branch_id, other in adjacency[node]:
            if other not in path:
                path[other] = path[node] + z_line.get(branch_id, 0.0)
                queue.append(other)
    assert path.keys() == network.node_u_rated.keys(), "a node is not reachable from the source"
    return path


def swer_of_nodes(network: Network) -> dict[int, Transformer1P]:
    """the ISO that feeds each SWER node (through one-phase lines from its secondary)"""
    adjacency = defaultdict(list)
    for line in network.lines:
        if line.phases == 1:
            adjacency[line.from_node].append(line.to_node)
            adjacency[line.to_node].append(line.from_node)
    iso_of: dict[int, Transformer1P] = {}
    for iso in sorted(network.isos, key=lambda t: t.id):
        iso_of[iso.to_node] = iso
        queue = deque([iso.to_node])
        while queue:
            node = queue.popleft()
            for other in adjacency[node]:
                if other not in iso_of:
                    iso_of[other] = iso
                    queue.append(other)
    swer_nodes = {node for node, node_class in network.node_class.items() if node_class == CLASS_SWER}
    assert iso_of.keys() == swer_nodes, "a SWER node is not fed by exactly one ISO"
    return iso_of


def largest(items, key):
    """the item with the largest key; on a tie the one with the lowest id (items are (id, value) pairs)"""
    return min(items, key=lambda item: (-key(item), item[0]))


def fault_scenarios(network: Network) -> list[FaultScenario]:
    """the fault scenarios of the short-circuit case (the table of the README), selected on the merged network, so every
    fault node exists in both tools. Path impedance: the sum |z1| L of the lines on the breadth-first tree path
    (path_impedances); every tie goes to the lowest id"""
    path = path_impedances(network)
    iso_of = swer_of_nodes(network)
    scenarios: list[FaultScenario] = []

    def add(fault: tuple[int, str, str, float, float], location: str, reason: str):
        """fault: (node, fault_type, fault_phase, r_f, x_f)"""
        node, fault_type, phase, r_f, x_f = fault
        scenario = ScenarioFault(FAULT_ID, node, fault_type, phase, r_f, x_f, scenario=len(scenarios))
        scenarios.append(FaultScenario(scenario, location, reason))

    zone_bus = network.zone.to_node
    add((zone_bus, "three_phase", "abc", 0.0, 0.0), "zone-substation 22 kV bus", "the secondary of the zone substation")

    mv_nodes = [(node, path[node]) for node, node_class in network.node_class.items() if node_class == CLASS_MV]
    node, z = largest(mv_nodes, key=lambda item: item[1])
    add(
        (node, "single_phase_to_ground", "a", 0.0, 0.0),
        "22 kV line end",
        f"the 22 kV node with the largest path impedance from the zone-substation bus ({z - path[zone_bus]:.4g} ohm)",
    )

    regulator = min(network.regulators, key=lambda t: (-t.sn, t.id))
    add(
        (regulator.to_node, "two_phase", "bc", TWO_PHASE_R_F, 0.0),
        f"secondary bus of the regulator {network.regulator_names[regulator.id]}",
        (
            f"the largest regulator ({regulator.sn / 1e3:g} kVA); r_f = {TWO_PHASE_R_F:g} ohm, since a bolted 22 kV"
            " bridge has a round-off of about 1.6e-2 in OpenDSS"
        ),
    )

    swer_nodes = [(node, path[node] - path[iso.to_node]) for node, iso in iso_of.items()]
    node, z = largest(swer_nodes, key=lambda item: item[1])
    iso_name = network.iso_names[iso_of[node].id]
    reason = f"the SWER node with the largest path impedance from its ISO secondary ({iso_name}, {z:.4g} ohm)"
    add((node, "single_phase_to_ground", "a", 0.0, 0.0), "SWER line end", reason)
    add(
        (node, "single_phase_to_ground", "a", SWER_FAULT_Z.real, SWER_FAULT_Z.imag),
        "SWER line end (high-impedance fault, a Reactor in OpenDSS)",
        "the node of scenario 3",
    )

    def by_kind(kind: str) -> list[Transformer1P | ThreePhaseTransformer]:
        return [t for t in network.distribution if network.distribution_kind[t.id] == kind]

    swer_transformers = [(t.id, path[t.from_node] - path[iso_of[t.from_node].to_node], t) for t in by_kind(KIND_SWER)]
    _, z, t = largest(swer_transformers, key=lambda item: item[1])
    add(
        (t.to_node, "single_phase_to_ground", "a", 0.0, 0.0),
        "LV busbar of a SWER transformer (240 V)",
        (
            f"the SWER transformer {t.id} ({network.distribution_kva[t.id]:g} kVA) whose primary has the largest path"
            f" impedance from its ISO secondary ({network.iso_names[iso_of[t.from_node].id]}, {z:.4g} ohm)"
        ),
    )

    for kind, fault_type, phase, location in (
        (KIND_SPLIT, "two_phase", "ab", "LV busbar of a split-phase transformer (250-0-250 V)"),
        (KIND_3PH, "three_phase", "abc", "LV busbar of a three-phase transformer (433 V)"),
        (KIND_1PH_22KV, "single_phase_to_ground", "a", "LV busbar of a 22 kV phase-to-earth transformer (250 V)"),
    ):
        candidates = [(t.id, network.distribution_kva[t.id], t) for t in by_kind(kind)]
        if not candidates:
            continue
        _, kva, t = largest(candidates, key=lambda item: item[1])
        reason = f"the largest {kind} transformer, {t.id} ({kva:g} kVA)"
        if kind == KIND_SPLIT:
            reason += "; a bolted two-phase fault is accurate in OpenDSS here (the shorted legs sit near 0 V)"
        add((t.to_node, fault_type, phase, 0.0, 0.0), location, reason)
    return scenarios


def three_phase_sc_reference(t: ThreePhaseTransformer) -> dict:
    """per-phase currents of both sides of a PGM transformer (the neutral conductors dropped), with the mirror
    Reactor's current on the from side"""
    i_from, i_to = terminal_values(f"Transformer.t{t.id}", dss.CktElement.Currents, 4)
    if t.low_susceptance() is not None:
        i_low, _ = terminal_values(f"Reactor.lo{t.id}n{t.from_node}", dss.CktElement.Currents, 3)
        i_from = i_from + i_low
    i_from_mag, i_from_angle = mag_angle(i_from)
    i_to_mag, i_to_angle = mag_angle(i_to)
    return {"id": t.id, "i_from": i_from_mag, "i_from_angle": i_from_angle, "i_to": i_to_mag, "i_to_angle": i_to_angle}


def sc_reference(network: Network, i_f: np.ndarray) -> dict[str, list[dict]]:
    """the reference of a solved short-circuit scenario: fault, single_phase_transformer and transformer, each sorted
    by id and rounded to 12 significant digits"""
    i_f_mag, i_f_angle = mag_angle(i_f)
    return rounded(
        {
            "fault": [{"id": FAULT_ID, "i_f": i_f_mag, "i_f_angle": i_f_angle}],
            "single_phase_transformer": [
                transformer_sc_reference(t, None)
                for t in sorted(network.single_phase_transformers(), key=lambda t: t.id)
            ],
            "transformer": [
                three_phase_sc_reference(t) for t in sorted(network.three_phase_transformers(), key=lambda t: t.id)
            ],
        }
    )


def dss_sc_setup(network: Network) -> list[str]:
    """the network with the IEC 60909 source (pu = c, impedances times c), no loads, and the solution settings"""
    return [
        *dss_network(network, voltage_pu=SC_VOLTAGE_FACTOR),
        f"Set VoltageBases={DSS_VOLTAGE_BASES}",
        "CalcVoltageBases",
        f"Set Tolerance={DSS_TOLERANCE:g}",
        f"Set MaxIterations={DSS_MAX_ITERATIONS}",
        "Set Mode=Snapshot",
    ]


def fault_definitions(fault: ScenarioFault, bridge_scale: float = 1.0) -> list[str]:
    return [f"New {element} {definition}" for element, definition in fault.dss_elements(bridge_scale)]


def solve_fault(network: Network, commands: list[str], fault: ScenarioFault) -> np.ndarray:
    """run the commands (ending in Solve) and return the fault current"""
    run_commands(commands)
    assert dss.Solution.Converged(), f"{network.spec.name} scenario {fault.scenario}: OpenDSS did not converge"
    return fault_currents(fault)


@dataclass
class ShortCircuitRun:
    script: list[str]  # test.dss: exactly the commands that were run
    scenarios: list[FaultScenario]
    references: list[dict]  # per scenario
    currents: list[np.ndarray]  # the fault current per scenario
    iterations: list[int]
    bridge: dict[int, float] = field(default_factory=dict)  # scenario -> the change of i_f with a 10x bridge
    reactor: dict[int, float] = field(default_factory=dict)  # scenario -> the deviation of the Reactor from V / z_f
    source_currents: list[np.ndarray] = field(default_factory=list)  # the fault currents in the source network


def run_short_circuit(network: Network, source: Network) -> ShortCircuitRun:
    """the batch short circuit, as test.dss runs it: the network once (IEC 60909 maximum, no loads), then per scenario
    the fault elements, Solve, the reference read and the fault elements disabled. Self-check 6: a Reactor fault
    carries V / z_f, and a bolted two-phase fault (a 1e-13 ohm bridge) gives the same fault current with a 10x larger
    bridge, solved on its own afterwards. Finally the same faults on the source network (every line, the sheet uk of
    the zone substation), for the effect of the merge"""
    name = network.spec.name
    scenarios = fault_scenarios(network)
    dss.Basic.ClearAll()
    run = ShortCircuitRun(
        script=[
            (
                f"! short circuit (IEC 60909 maximum, c = {SC_VOLTAGE_FACTOR:g}) on the Australian SWER network {name};"
                " written by tests/data/generators/australian_swer_networks.py, which runs these commands line by line"
            ),
            "! per scenario: the fault elements, Solve, (the reference is read), the fault elements disabled",
        ],
        scenarios=scenarios,
        references=[],
        currents=[],
        iterations=[],
    )
    run_commands(run.script)
    setup = dss_sc_setup(network)
    run_commands(setup)
    run.script += setup
    for k, scenario in enumerate(scenarios):
        fault = scenario.fault
        commands = [
            (
                f"! scenario {k}: {fault.fault_type} {fault.fault_phase} at node {fault.node} ({scenario.location}),"
                f" r_f = {fault.r_f:g} ohm, x_f = {fault.x_f:g} ohm"
            ),
            *fault_definitions(fault),
            "Solve",
        ]
        i_f = solve_fault(network, commands, fault)
        run.iterations.append(dss.Solution.Iterations())
        if fault.x_f != 0.0:
            check_fault_impedance(fault, i_f)
            assert fault.fault_type == "single_phase_to_ground", "the README value covers one faulted phase"
            k_phase = "abc".index(fault.fault_phase)
            expected = bus_phasors(fault.node)[k_phase + 1] / complex(fault.r_f, fault.x_f)
            run.reactor[k] = float(abs(i_f[k_phase] - expected) / abs(expected))
        run.currents.append(i_f)
        run.references.append(sc_reference(network, i_f))
        disable = [f"Edit {element} enabled=no" for element, _ in fault.dss_elements()]
        run_commands(disable)
        run.script += [*commands, *disable]
        print(
            f"  scenario {k}: {fault.fault_type} {fault.fault_phase} at node {fault.node} ({scenario.location}):"
            f" |i_f| = {', '.join(f'{abs(v):.6g}' for v in i_f)} A ({run.iterations[-1]} iterations)"
        )
    for k, scenario in enumerate(scenarios):
        fault = scenario.fault
        if not fault.has_bridge():
            continue
        dss.Basic.ClearAll()
        i_f = solve_fault(network, [*setup, *fault_definitions(fault, bridge_scale=10.0), "Solve"], fault)
        dev = float(np.max(np.abs(i_f - run.currents[k])) / np.max(np.abs(run.currents[k])))
        assert dev < 1e-10, f"{name} scenario {k}: the bolted two-phase fault is not accurate ({dev:.1e} relative)"
        run.bridge[k] = dev
        print(f"  scenario {k}: a 10x larger bridge changes i_f by {dev:.1e} relative")
    dss.Basic.ClearAll()
    run_commands(dss_sc_setup(source))
    for k, scenario in enumerate(scenarios):
        fault = scenario.fault
        i_f = solve_fault(source, [*fault_definitions(fault), "Solve"], fault)
        run_commands([f"Edit {element} enabled=no" for element, _ in fault.dss_elements()])
        run.source_currents.append(i_f)
        print(f"  scenario {k}: against the source network, i_f changes by {source_change(run, k):.2e} relative")
    return run


def source_change(run: ShortCircuitRun, k: int) -> float:
    """the largest change of a fault current phasor (relative to the largest one) of the merged network against the
    source network"""
    reference = run.source_currents[k]
    return float(np.max(np.abs(run.currents[k] - reference)) / np.max(np.abs(reference)))


# ---------------------------------------------------------------------------------------------------------------
# writing


def compact_dataset(dataset_type: str, data, *, is_batch: bool) -> str:
    """a PGM dataset in the compact JSON format: the root "attributes" lists the attributes of each component (the
    keys of its rows in the order of their first appearance; a row without one of them gets null, the absent value)
    and each row is a list, written one per line"""
    scenarios = data if is_batch else [data]
    attributes: dict[str, list[str]] = {}
    for scenario in scenarios:
        for component, rows in scenario.items():
            names = attributes.setdefault(component, [])
            for row in rows:
                names += [name for name in row if name not in names]
    attributes = {component: names for component, names in attributes.items() if names}

    def scenario_text(scenario) -> str:
        blocks = []
        for component, rows in scenario.items():
            if not rows:
                continue
            body = ",\n".join(
                json.dumps([row.get(name) for name in attributes[component]], separators=(",", ":")) for row in rows
            )
            blocks.append(f"{json.dumps(component)}:[\n{body}\n]")
        return "{\n" + ",\n".join(blocks) + "\n}"

    header = json.dumps(
        {"version": "1.0", "type": dataset_type, "is_batch": is_batch, "attributes": attributes},
        separators=(",", ":"),
    )
    if is_batch:
        body = "[\n" + ",\n".join(scenario_text(scenario) for scenario in scenarios) + "\n]"
    else:
        body = scenario_text(data)
    return f'{header[:-1]},\n"data":{body}}}\n'


def write_licensed(path: Path, text: str):
    path.write_text(text, encoding="utf-8")
    Path(f"{path}.license").write_text(LICENSE, encoding="utf-8")


def load_update(loads: list[Load]) -> dict[str, list[dict]]:
    return {
        "asym_load": [
            {"id": load.id, "p_specified": list(load.p), "q_specified": list(load.q)}
            for load in sorted(loads, key=lambda load: load.id)
        ]
    }


def count_table(network: Network, with_loads: bool) -> list[str]:
    mv = sum(1 for node_class in network.node_class.values() if node_class == CLASS_MV)
    swer = sum(1 for node_class in network.node_class.values() if node_class == CLASS_SWER)
    n_three = sum(1 for line in network.lines if line.phases == 3)
    kinds = Counter(network.distribution_kind.values())
    regulators = ", ".join(
        f"{network.regulator_names[t.id]} ({t.sn / 1e3:g} kVA)" for t in sorted(network.regulators, key=lambda t: t.id)
    )
    isos = ", ".join(
        f"{network.iso_names[t.id]} ({t.sn / 1e3:g} kVA)" for t in sorted(network.isos, key=lambda t: t.id)
    )
    rows = [
        (
            "nodes",
            (
                f"{len(network.node_u_rated)}: the 66 kV source bus, {mv} 22 kV, {swer} SWER and"
                f" {len(network.distribution)} LV busbars"
            ),
        ),
        ("`line` (three-phase, 22 kV)", f"{n_three}"),
        ("`asym_line` (one-phase SWER, phase a)", f"{len(network.lines) - n_three}"),
        ("zone substation (`transformer`)", f"66/{network.zone.u2 / 1e3:g} kV, {network.zone.sn / 1e6:g} MVA, Dyn1"),
        ("voltage regulators (`transformer`, YNyn0)", f"{len(network.regulators)}: {regulators}"),
        ("ISO transformers (`single_phase_transformer`)", f"{len(network.isos)}: {isos}"),
        ("three-phase distribution transformers (`transformer`, 22/0.433 kV Dyn1)", f"{kinds[KIND_3PH]}"),
        (
            "split-phase distribution transformers (`single_phase_transformer`, 22 kV / 250-0-250 V)",
            f"{kinds[KIND_SPLIT]}",
        ),
        ("SWER distribution transformers (`single_phase_transformer`, 12.7/0.24 kV)", f"{kinds[KIND_SWER]}"),
        (
            "22 kV phase-to-earth distribution transformers (`single_phase_transformer`, 12.7/0.25 kV)",
            f"{kinds[KIND_1PH_22KV]}",
        ),
        ("capacitors (`shunt`, 900 kvar)", f"{len(network.capacitors)}"),
    ]
    if with_loads:
        rows.append(("loads (`asym_load`, const_pq)", f"{len(network.loads_min)}"))
    return [
        "| Element | Count |",
        "|---|---|",
        *(f"| {element} | {count} |" for element, count in rows if count != "0"),
    ]


def results_table(network: Network, summaries: list[PowerFlowSummary]) -> list[str]:
    header = "| | " + " | ".join(f"scenario {k}: {s.label}" for k, s in enumerate(summaries)) + " |"
    lines = [header, "|---" * (len(summaries) + 1) + "|"]

    def row(name, values):
        lines.append(f"| {name} | " + " | ".join(values) + " |")

    row("OpenDSS iterations", [str(s.iterations) for s in summaries])
    row("feeder P, Q", [f"{s.feeder_p / 1e6:.3f} MW, {s.feeder_q / 1e6:.3f} Mvar" for s in summaries])
    row("load P", [f"{s.load_p / 1e6:.3f} MW" for s in summaries])
    row("losses", [f"{s.losses / 1e3:.1f} kW" for s in summaries])
    for node_class, nominal in NOMINAL_VOLTAGE.items():
        if node_class in summaries[0].extremes:
            row(
                f"{node_class} (pu of {nominal:.1f} V)",
                [f"{s.extremes[node_class][0]:.4f} to {s.extremes[node_class][2]:.4f}" for s in summaries],
            )
    row(
        "highest ISO loading",
        [
            "{1:.1%} ({0})".format(
                *max(((network.iso_names[i], v) for i, v in s.iso_loading.items()), key=lambda x: x[1])
            )
            for s in summaries
        ],
    )
    row(
        "highest regulator loading",
        [
            "{1:.1%} ({0})".format(
                *max(((network.regulator_names[i], v) for i, v in s.regulator_loading.items()), key=lambda x: x[1])
            )
            for s in summaries
        ],
    )
    return lines


PF_README_HEAD = """\
<!--
{license}-->

# Validation case: power flow on the Australian rural SWER network {title}

Asymmetric power flow ({methods}) on the medium-voltage part of the rural 22 kV feeder {title} of AusNet Services
(Victoria, Australia), with its SWER (single-wire earth return) networks, with OpenDSS as the reference solver. It is a
batch of two scenarios: scenario 0 at minimum and scenario 1 at maximum demand. `input.json` holds the minimum demand
and `update_batch.json` sets the `asym_load` powers of both scenarios. There is no single-calculation output
(`asym_output.json`).

"""

# the parts shared by the power-flow and the short-circuit README
README_SHARED = """\
## Source data

- Repository: {repo}, commit `{commit}`
- File: `{file}`, sha256 `{sha256}`
- Licence: BSD 3-Clause, {copyright} (`LICENSES/BSD-3-Clause.txt`). The files of this case are
  derived from it and are licensed `MPL-2.0 AND BSD-3-Clause` (see the `.license` files).
- The model is derived from the repository's `backend.py` (its OpenDSS model of the network), with the deviations
  listed below. The workbook is not part of this repository.

## Network

{counts}

## Modelling and deviations from backend.py

- The LV networks are replaced by one lumped load per distribution transformer on a new LV busbar. The medium-voltage
  network, the SWER networks and every transformer of the source data are kept.
- Voltage regulators (every `mvtx` row with equal primary and secondary kV, as in backend; on KLO14 this includes
  PLENTY_VLY_RISE, 22/22 kV and 200 kVA) are fixed at the neutral tap as a YNyn0 bank on the rating with backend's xhl
  and %loadloss. backend builds each one as three 12.7/1.27 kV autotransformers under a RegControl starting at tap 0;
  at the neutral tap they are a 1:1 connection plus the leakage impedance, which both tools model identically as
  YNyn0 (in positive and zero sequence).
- ISO transformers (22/12.7 kV) are two-winding single-phase transformers from the phase pair of `Conn_Type` (R, W, B =
  a, b, c) to phase a of the SWER network, with clock 6 if `Conn_Type` lists the pair in reverse order; backend's open
  third winding is dropped. The OpenDSS buses are written literally from `Conn_Type`, so OpenDSS independently checks
  the clock-6 mapping.
- Split-phase distribution transformers (22 kV phase-to-phase, 250-0-250 V) are centre-tapped single-phase
  transformers (u2 = 500 V) with the Kersting default split of uk and pk; OpenDSS gets the matching XHL = XHT and XLT
  instead of backend's default XHT = 35 % and XLT = 30 %.
- Three-phase distribution transformers (22/0.433 kV) are Dyn1, OpenDSS's default delta-wye (the `Conn_Type` rotation
  is ignored, as in backend).
- KLO14 has 44 distribution transformers with a 22/0.433 kV rating that are wired from phase R to earth. They are
  12.7/0.25 kV single-phase transformers here (both ratings taken phase-to-earth, which keeps the data's turns ratio).
- Every transformer has i0 = p0 = 0 (the no-load losses and the magnetising current are dropped): PGM places the
  magnetising branch as a pi, OpenDSS at winding 1. uk and pk come from the sheet's xhl and loadloss, and OpenDSS gets
  XHL and %loadloss computed back from them.
- Distribution transformer taps are written explicitly (tap_pos = (wdg1_tap - 1) / 0.025 on the from side, tap_size =
  u1 / 40, limits -2 to 2). backend's maxtap = 1.137 and mintap = 1.028 are ignored, as in backend, where the tap stays
  at wdg1_tap.
- PGM adds a low susceptance (1e-8 sn / uk) to the zero sequence of a delta winding and to the common mode of a
  phase-to-phase single-phase winding. OpenDSS gets matching shunt Reactors (`Reactor.lo{{id}}n{{node}}`), whose powers
  are added to the from side of the transformer in the reference.
- The source has Z1 = 0.52824 + j2.113 ohm and Z0 with the magnitude of backend's 0.59157 + j1.7 ohm but Z1's angle
  (PGM has a real `z01_ratio`). This is irrelevant: the zone substation primary is delta, so the 66 kV source carries
  no zero-sequence current.
- Line rows named `Delete` are dropped. One-phase SWER lines use Z1 and C1 of their LineCode only, as OpenDSS solves
  backend's 1-phase LineCodes (r0, x0 and b0 are ignored), so they have no earth-return reactance in either tool. All
  22 kV and 12.7 kV SWER nodes have `u_rated` 22 kV (the SWER convention of power-grid-model).
- Capacitors are grounded-wye shunts with b1 = b0 = kvar / kV^2.
- Near-zero lines are merged (in both tools): see below.
- OpenDSS runs with tolerance {dss_tolerance} (1e-12 stalls at round-off on SMR8) and MaxIterations
  {dss_max_iterations}.

## Merged near-zero lines

The source data has connector lines with a series impedance down to 1.4e-7 ohm and a zone substation with uk =
1.41e-6. They make the admittance matrix so ill-conditioned that PGM's Newton-Raphson stalls at a max deviation of
about 3e-7 (error tolerance 1e-8) and PGM and OpenDSS agree to only about 2e-7. This is not caused by the reduced-phase
or SWER parts: the three-phase part of the network on its own (the upstream code path) stalls at the same level.

Therefore every line whose series impedance |z1| L is below {z_merge} ohm ({z_merge_pu} pu on the 22 kV base of 484 ohm,
1 MVA) is removed and its two end nodes are merged into one, in both tools. A merged group of nodes keeps the node
closest to the source (the fewest branches; on a tie the lowest id), and every element on a merged node moves to it:
lines, transformer terminals and capacitors. A three-phase line merges 22 kV nodes and a one-phase line SWER nodes;
the generator asserts that a 22 kV node never merges with a SWER node. A kept line whose ends merge is dropped (it runs
in parallel with merged lines only), and lines that end up between the same two nodes stay parallel. Every other line
keeps the length and LineCode of the data. The zone substation uk is raised to {zone_uk} (uk and pk scaled together,
r / x kept): with a smaller uk it is the smallest impedance left after the merge, and PGM misses the tolerances on
KLO14. Here: {merge}
"""

PF_README_TAIL = """
The merge changes the OpenDSS solution against the source network (every line of the data and the sheet uk of the zone
substation, solved with the same tolerance) by the following. A node voltage change is the largest
|V - V_source| / |V_source| over the phases of the kept nodes of a class (the phasors, so the angle counts too); a
merged node is compared at the node it is merged into.

{deviation}

## Demand

One const_pq `asym_load` per distribution transformer at power factor {power_factor} lagging, a fraction of its kVA:
{min_demand} at minimum and {max_demand} at maximum demand. Three-phase loads are split equally over the phases,
split-phase loads over the two legs (to the centre tap), and phase-to-earth loads are on phase a. At maximum demand the
SWER loads behind an ISO get the coincidence factor min(1, S_iso / ({max_demand} x the kVA of the SWER transformers
behind it)); applied here: {factors}. A higher maximum (the 0.60 of the original plan) does not converge in OpenDSS:
the voltage collapses on both networks.

## OpenDSS results

{results}

## Reference output

`asym_output_batch.json` has, per scenario, `node` (u, u_angle; 0 on absent phases), `asym_line` (p and q, from and
to), `single_phase_transformer` and `transformer` (p and q, from and to, and the loading), rounded to 12 significant
digits. It has no three-phase `line` output: the node voltages at both ends pin every line, and it keeps the data small.
All datasets are written in the compact JSON format (a root `attributes` list per component and list rows).

Tolerances: rtol {rtol} and atol {atol}, the smallest round values with a 10x margin over both scenarios, with the
atol of the power fields at most 1e-3 (W, var). The power fields have their own atol: a branch power is computed from
the voltages at its ends, so a short branch amplifies their round-off. {agreement}

## IDs

An MV node keeps its AusNet Node_ID; the 66 kV source bus is 1000000; the LV busbar of distribution transformer E is
2000000 + E; a line, transformer, capacitor or the zone substation (Element_ID E) is 3000000 + E; the load of
distribution transformer E is 4000000 + E; the source is 5000000. In OpenDSS every bus is `n{{node id}}`, and the
elements are `Line.l{{id}}`, `Transformer.t{{id}}`, `Capacitor.c{{id}}` and `Load.ld{{id}}_{{a,b,c}}`.

## Generation

The files are generated by `tests/data/generators/australian_swer_networks.py` with OpenDSS (opendssdirect.py
{opendss_version}) as the reference solver; do not edit them by hand. `test.dss` is exactly the script that was run:
the generator runs it line by line and reads each reference right after the `Solve` of its scenario. From the
repository root (the workbooks are downloaded from the pinned commit, or read from a local directory with
`--source-dir`):

```sh
uv run --no-project --with opendssdirect.py --with numpy --with openpyxl \\
    python tests/data/generators/australian_swer_networks.py tests/data {name}
```
"""

README_TEMPLATE = PF_README_HEAD + README_SHARED + PF_README_TAIL

SC_README_HEAD = """\
<!--
{license}-->

# Validation case: short circuit on the Australian rural SWER network {title}

IEC 60909 short circuit (maximum, c = {c}) on the medium-voltage part of the rural 22 kV feeder {title} of AusNet
Services (Victoria, Australia), with its SWER (single-wire earth return) networks, with OpenDSS as the reference
solver. The network is the one of the power-flow case `tests/data/power_flow/australian-swer/{name}`, without loads.
It is a batch of {n_scenarios} fault scenarios: `input.json` has the fault (one ID, {fault_id}) of scenario 0, and
`update_batch.json` points it to the node of each scenario (fault_object, fault_type, fault_phase, r_f, x_f). There is
no single-calculation output (`sc_output.json`).

"""

SC_README_TAIL = """
The faults are selected on the merged network, so every fault node exists in both tools (a merged node does not). The
merge changes the OpenDSS fault currents against the source network (every line of the data and the sheet uk of the
zone substation, with the same faults) by the following: the largest |I_f - I_f,source| over the phases, relative to
the largest |I_f,source|.

{source_change}

## Fault scenarios

One fault (ID {fault_id}) per scenario, selected deterministically by the generator. The path impedance of a node is
the sum |z1| L of the lines on its path from the source bus along the breadth-first tree of the network (the branches
of a node taken in the order of their IDs; transformers count 0); every tie goes to the lowest ID. A bolted fault
(r_f = 0) is a Fault of 1e-13 ohm in OpenDSS (OpenDSS replaces r = 0 by 1e-4 ohm), a fault with x_f a series R + jX
Reactor.

{faults}

Self-checks of the generator: OpenDSS converges in every scenario; the Reactor of the high-impedance fault carries
V / z_f with the OpenDSS node voltage ({reactor}); the bolted two-phase fault, solved again on its own with a 10x
larger bridge, gives the same fault current ({bridge}).

## Short-circuit model

- The source is c times 66 kV behind c times its impedance (IEC 60909, as in PGM), with c = {c} for the maximum
  (PGM takes it from the u_rated of the source node).
- No loads. The capacitors and the line capacitances stay in, as in PGM's short-circuit admittance matrix.
- `test.dss` builds the network once; per scenario it defines the fault elements (`Fault.f{fault_id}s{{k}}` or
  `Reactor.f{fault_id}s{{k}}` for scenario k), solves, and disables them again. OpenDSS runs with tolerance
  {dss_tolerance} in snapshot mode.

## Reference output

`sc_output_batch.json` has, per scenario, `fault` (i_f, i_f_angle), `single_phase_transformer` and `transformer` (i
and i_angle, from and to; the current of a low-susceptance Reactor is added to the from side), rounded to 12
significant digits. It has no `node`, `line` or `asym_line` output: {size_reason}. All datasets are written in the
compact JSON format (a root `attributes` list per component and list rows).

Tolerances: rtol {rtol} and atol {atol}, the smallest round values with a 10x margin over all scenarios. {agreement}

## IDs

As in the power-flow case: an MV node keeps its AusNet Node_ID; the 66 kV source bus is 1000000; the LV busbar of
distribution transformer E is 2000000 + E; a line, transformer, capacitor or the zone substation (Element_ID E) is
3000000 + E; the source is 5000000 and the fault {fault_id}. In OpenDSS every bus is `n{{node id}}`, and the elements
are `Line.l{{id}}`, `Transformer.t{{id}}` and `Capacitor.c{{id}}`.

## Generation

The files are generated by `tests/data/generators/australian_swer_networks.py` with OpenDSS (opendssdirect.py
{opendss_version}) as the reference solver, together with the power-flow case; do not edit them by hand. `test.dss` is
exactly the script that was run: the generator runs it line by line and reads each reference right after the `Solve`
of its scenario. From the repository root (the workbooks are downloaded from the pinned commit, or read from a local
directory with `--source-dir`):

```sh
uv run --no-project --with opendssdirect.py --with numpy --with openpyxl \\
    python tests/data/generators/australian_swer_networks.py tests/data {name}
```
"""

SC_README_TEMPLATE = SC_README_HEAD + README_SHARED + SC_README_TAIL
# measured: with asym_line output the four cases would take 10.0 MB, above the budget of 10 MB (8.97 MB without it)
SC_SIZE_REASON = (
    "the transformer currents pin the network around each fault, and the data of the four cases on both networks"
    " has a size budget of 10 MB (it is 9.0 MB; the asym_line output would add 1.05 MB and the node output 5 MB)"
)


def atol_text(atol) -> str:
    if not isinstance(atol, dict):
        return f"{atol:g}"
    others = ", ".join(f"{value:g} on `{pattern}`" for pattern, value in atol.items() if pattern != "default")
    return f"{atol['default']:g} ({others})"


def shared_readme_fields(network: Network, with_loads: bool) -> dict:
    """the fields of README_SHARED"""
    spec = network.spec
    return {
        "license": LICENSE,
        "title": spec.name.upper(),
        "name": spec.name,
        "repo": SOURCE_REPO,
        "commit": SOURCE_COMMIT,
        "file": spec.file,
        "sha256": spec.sha256,
        "copyright": UPSTREAM_COPYRIGHT,
        "counts": "\n".join(count_table(network, with_loads)),
        "dss_tolerance": f"{DSS_TOLERANCE:g}",
        "dss_max_iterations": DSS_MAX_ITERATIONS,
        "z_merge": f"{Z_MERGE_OHM:g}",
        "z_merge_pu": f"{Z_MERGE_OHM / (U_MV**2 / 1e6):.2g}",
        "zone_uk": f"{ZONE_SUBSTATION_UK:g}",
        "merge": merge_text(network),
        "opendss_version": dss.__version__,
    }


def readme_text(network: Network, run: PowerFlowRun, deviation: SourceDeviation, params: dict) -> str:
    methods = params["calculation_method"]
    methods = methods if isinstance(methods, list) else [methods]
    factors = [
        f"{network.iso_names[iso.id]} {network.iso_factor[iso.id]:.4g} ({iso.sn / 1e3:g} kVA for"
        f" {network.iso_swer_kva[iso.id]:g} kVA of SWER transformers)"
        for iso in sorted(network.isos, key=lambda t: t.id)
        if network.iso_factor[iso.id] < 1.0
    ]
    return README_TEMPLATE.format(
        **shared_readme_fields(network, with_loads=True),
        methods=", ".join(f"`{m}`" for m in methods),
        deviation="\n".join(deviation_table(deviation)),
        power_factor=POWER_FACTOR,
        min_demand=f"{MIN_DEMAND:g}",
        max_demand=f"{MAX_DEMAND:g}",
        factors="; ".join(factors) if factors else "none (no ISO is above its rating)",
        results="\n".join(results_table(network, run.summaries)),
        rtol=f"{params['rtol']:g}",
        atol=atol_text(params["atol"]),
        agreement=PF_AGREEMENT[network.spec.name],
    )


def fault_table(run: ShortCircuitRun) -> list[str]:
    lines = [
        "| # | node | where | fault | r_f, x_f (ohm) | why this node | OpenDSS \\|I_f\\| (A) |",
        "|---|---|---|---|---|---|---|",
    ]
    for k, (scenario, i_f) in enumerate(zip(run.scenarios, run.currents)):
        fault = scenario.fault
        magnitudes = ", ".join(f"{abs(i_f['abc'.index(p)]):.6g}" for p in fault.fault_phase)
        lines.append(
            f"| {k} | {fault.node} | {scenario.location} | {fault.fault_type} {fault.fault_phase} |"
            f" {fault.r_f:g}, {fault.x_f:g} | {scenario.reason} | {magnitudes} |"
        )
    return lines


def source_change_table(run: ShortCircuitRun) -> list[str]:
    header = "| scenario | " + " | ".join(str(k) for k in range(len(run.scenarios))) + " |"
    return [
        header,
        "|---" * (len(run.scenarios) + 1) + "|",
        "| change of I_f | " + " | ".join(f"{source_change(run, k):.1e}" for k in range(len(run.scenarios))) + " |",
    ]


def sc_readme_text(network: Network, run: ShortCircuitRun, params: dict) -> str:
    reactor = "; ".join(f"scenario {k}: {dev:.1e} relative" for k, dev in run.reactor.items())
    bridge = "; ".join(f"scenario {k}: {dev:.1e} relative" for k, dev in run.bridge.items())
    return SC_README_TEMPLATE.format(
        **shared_readme_fields(network, with_loads=False),
        c=f"{SC_VOLTAGE_FACTOR:g}",
        n_scenarios=len(run.scenarios),
        fault_id=FAULT_ID,
        source_change="\n".join(source_change_table(run)),
        faults="\n".join(fault_table(run)),
        reactor=reactor,
        bridge=bridge,
        size_reason=SC_SIZE_REASON,
        rtol=f"{params['rtol']:g}",
        atol=atol_text(params["atol"]),
        agreement=SC_AGREEMENT[network.spec.name],
    )


def write_power_flow_case(network: Network, run: PowerFlowRun, deviation: SourceDeviation, root: Path):
    target = root / "power_flow" / "australian-swer" / network.spec.name
    target.mkdir(parents=True, exist_ok=True)
    params = PF_PARAMS[network.spec.name]
    write_licensed(
        target / "input.json", compact_dataset("input", pgm_rows(network, network.loads_min), is_batch=False)
    )
    write_licensed(
        target / "update_batch.json",
        compact_dataset("update", [load_update(network.loads_min), load_update(network.loads_max)], is_batch=True),
    )
    write_licensed(target / "asym_output_batch.json", compact_dataset("asym_output", run.references, is_batch=True))
    write_licensed(target / "params.json", json.dumps(params, indent=2) + "\n")
    write_licensed(target / "test.dss", "\n".join(run.script) + "\n")
    (target / "README.md").write_text(readme_text(network, run, deviation, params), encoding="utf-8")
    print(f"  wrote {target}")


def write_short_circuit_case(network: Network, run: ShortCircuitRun, root: Path):
    target = root / "short_circuit" / "australian-swer" / network.spec.name
    target.mkdir(parents=True, exist_ok=True)
    params = SC_CASE_PARAMS[network.spec.name]
    rows = pgm_rows(network, [])
    rows["fault"] = [pgm_fault(run.scenarios[0].fault, with_status=True)]
    write_licensed(target / "input.json", compact_dataset("input", rows, is_batch=False))
    write_licensed(
        target / "update_batch.json",
        compact_dataset(
            "update", [{"fault": [pgm_fault(s.fault, with_status=False)]} for s in run.scenarios], is_batch=True
        ),
    )
    write_licensed(target / "sc_output_batch.json", compact_dataset("sc_output", run.references, is_batch=True))
    write_licensed(target / "params.json", json.dumps(params, indent=2) + "\n")
    write_licensed(target / "test.dss", "\n".join(run.script) + "\n")
    (target / "README.md").write_text(sc_readme_text(network, run, params), encoding="utf-8")
    print(f"  wrote {target}")


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
    specs = [spec for spec in NETWORKS if not args.networks or spec.name in args.networks]
    print(f"opendssdirect.py {dss.__version__}; MAX_DEMAND = {MAX_DEMAND:g}, MIN_DEMAND = {MIN_DEMAND:g}")
    with tempfile.TemporaryDirectory() as download_dir:
        paths = fetch_sources(specs, args.source_dir, Path(download_dir))
        sheets = {spec.name: read_workbook(paths[spec.name]) for spec in specs}
    for spec in specs:
        source = build_network(spec, sheets[spec.name])
        network = merge_connectors(source)
        check_elements(network)
        run = run_power_flow(network)
        deviation = source_deviation(network, source)
        sc_run = run_short_circuit(network, source)
        if not args.check_only:
            write_power_flow_case(network, run, deviation, args.root)
            write_short_circuit_case(network, sc_run, args.root)


if __name__ == "__main__":
    main()
