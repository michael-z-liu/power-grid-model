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
  missing-phase-for-winding, split-phase, split-phase-flat-start, split-phase-linear}
- tests/data/short_circuit/{single-phase-transformer-swer, single-phase-transformer-split-phase}
- tests/data/short_circuit/{swer-fault-impedance, swer-multiple-faults, swer-c-minimum, swer-fault-on-absent-phase}
- tests/data/power_flow/asym-load-delta/{three-phase, three-phase-linear, two-phase, two-phase-linear, split-phase,
  split-phase-flat-start, split-phase-linear, symmetric-not-supported, delta-on-absent-phase}
- tests/data/state_estimation/asym-load-delta-not-supported

Run from the repository root (opendssdirect.py is not a project dependency, so it is added on the fly):

    uv run --no-project --with opendssdirect.py --with numpy \\
        python tests/data/generators/reduced_phase_opendss_cases.py tests/data

Optionally pass case path fragments after the root to generate only the matching cases.

Modelling choices:

- source: Vsource with Z1/Z0 in ohm from sk, rx_ratio and z01_ratio, exactly as PGM computes them; in a short
  circuit the voltage is c * u_rated (pu=c) and the impedance is c * u_rated^2 / sk (IEC 60909, as in PGM). c is
  taken as PGM takes it, from the u_rated of the source node (not of the fault location): 1.1 for maximum voltage
  scaling; for minimum 0.95 if the source node is at most 1000 V and 1.0 otherwise
- faults: OpenDSS replaces a fault resistance of 0 by 1e-4 ohm, so a bolted fault (r_f = 0) gets 1e-13 ohm instead,
  which is below the OpenDSS agreement of about 3e-10 relative. A bolted two-phase fault is a 1e-13 ohm bridge
  between the phases; its 1e13 S conductance swamps the nodal admittance of the two phases in double precision and
  leaves a round-off error of about eps |V| G in the fault current. This is harmless only where the shorted phases sit
  near 0 V, as at a split-phase node, because then eps |V| G is negligible; it is not harmless at a low-impedance node
  in general (1.3e-4 relative at a three-phase 400 V node, 1.6e-2 at 22 kV). Each bolted two-phase fault is therefore
  solved again with a 10x larger bridge and must agree to 1e-10 relative, which catches it. A bolted
  two-phase-to-ground fault (PGM: both phases at 0 V) is two Fault elements of 1e-13 ohm from each phase to earth,
  without a bridge. With r_f > 0, PGM joins the two phases and connects them to earth through z_f; in OpenDSS this
  needs a bridge between the phases at a node with a non-zero voltage, where the same round-off limits the reference
  to about 3e-9 relative (LV) and 1e-8 relative (22 kV) for any bridge resistance, too coarse for the tolerances of
  these cases, so it is not generated. The fault current per phase is the sum over all OpenDSS elements of the fault
  and over their terminals of the current from the node into the element
- fault impedance: an OpenDSS Fault is resistive only, so a fault with x_f != 0 (x_f > 0 only; OpenDSS may treat a
  Reactor with X = 0 specially) is a series R + jX Reactor from the faulted phase to earth (bus2 = n.0), between the
  two phases of a two-phase fault, or from each phase to earth for a three-phase fault. After the solve its current
  must equal V / (r_f + j x_f) with the OpenDSS node voltages to 1e-9 relative (check_fault_impedance, the decisive
  check of the Reactor R/X convention)
- several faults per scenario: an entry of Case.faults is a Fault or a tuple of simultaneous Faults. The input has
  every fault ID (the faults of the first scenario with status 1, the others with status 0 at their first
  appearance); with more than one ID the update sets each scenario's faults to status 1 and switches the others off,
  and the reference has i_f = 0 for a fault that is off (PGM's null output)
- lines: explicit R/X/C matrices of the present phases; a line with a neutral gets the numpy Kron-reduced matrix,
  which is cross-checked against an OpenDSS LineCode with Kron=yes; the exception is a line with opendss_kron=True,
  which gets the full matrices including the neutral as a LineCode with Kron=yes, so OpenDSS reduces the neutral itself
- loads: every PGM load (also a sym_load) becomes single-phase wye loads with kV = u_rated / sqrt(3) exactly;
  const_pq -> model=1, const_y -> model=2, const_i -> model=5; a generator is a load with negative power
- delta loads (connection = delta, M7): a 1-phase Load with conn=delta between the two phases of each non-zero entry
  (bus1=n.1.2, n.2.3 and n.3.1 for ab, bc and ca), or one native 3-phase conn=delta Load (bus1=n.1.2.3, kW and kvar
  the total of three equal entries); kV is the rated branch voltage, which PGM uses as U_delta,nom (decision D29):
  u_rated, or 2 u_rated / sqrt(3) (480 V) at a split-phase node (the to nodes of centre-tapped windings, spread along
  lines without phase c, as PGM does). After the solve the total power of every delta Load element must equal
  sum_k S_k (|V_p - V_q| / kV)^n to 1e-9 relative (check_delta_load_bases, the decisive check of the kV base and of the
  scaling of model 5 and 2; the 3-phase split is proven end-to-end by the per-phase load outputs). In the delta
  cases (load_output) the reference also has the
  asym_load and asym_gen outputs: p, q and i per phase are the sums per node of the terminal powers and currents of
  the Load elements of a load (negated for a generator), which are the phase quantities PGM reports (decision D30)
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
- centre-tapped (split-phase) transformers: a 1-phase, 3-winding wye-wye-wye Transformer (the primary and the two
  halves; half 1 from node 1 to the centre tap and half 2 from the centre tap to node 2, both reversed for clock 6)
  with kVs = u1, u2 / 2, u2 / 2, equal kVAs = sn, the star resistances of the primary and half branches z_0 and z_b
  per winding (%Rs; %loadloss is not used with three windings) and XHL = XHT = Im z_hs, XLT = 2 Im z_b on sn, where
  z_hs comes from uk_half/pk_half on sn / 2 or the Kersting default 1.5 r + j1.2 x. The centre tap is node 0, or node
  4 with an electrode Reactor to earth. A centre-tapped side has no low susceptance (a phase-to-phase primary keeps
  it). Each such transformer is first checked on its own: its OpenDSS YPrim must equal the numpy N-winding model
  A^T Z_B^-1 A (the decisive check of the three-winding conventions), and in open circuit U_b = -U_a with
  |U_a| = u2 / 2. The to-side powers are the sums per node over both halves, the loading is
  max(|sum s_from|, 2 |s_a|, 2 |s_b|) / sn (each half is rated sn / 2), and the centre-tap electrode current is checked
  to equal I_a + I_b (KCL at node 4)
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


# the phases (OpenDSS nodes) of the delta entries ab, bc and ca
DELTA_BRANCHES = ((1, 2), (2, 3), (3, 1))
DELTA_NAMES = ("ab", "bc", "ca")


@dataclass
class Load:
    id: int
    node: int
    kind: str  # asym_load, asym_gen, sym_load
    type: str
    p: list  # per phase (W) for asym, total for sym; for a delta connection the ab, bc and ca entries
    q: list
    status: int = 1
    connection: str = "wye"  # "wye" or "delta" (asym_load and asym_gen only)
    # a delta load with three equal entries as one native 3-phase OpenDSS delta load (instead of 1-phase loads)
    dss_three_phase: bool = False

    def delta(self):
        return self.connection == "delta"

    def pgm(self):
        result = {"id": self.id, "node": self.node, "status": self.status, "type": LOAD_TYPE[self.type]}
        if self.kind == "sym_load":
            result["p_specified"] = self.p
            result["q_specified"] = self.q
        else:
            result["p_specified"] = list(self.p)
            result["q_specified"] = list(self.q)
        if self.delta():
            assert self.kind != "sym_load", "a delta connection is only supported for asym_load and asym_gen"
            result["connection"] = 1
        return result

    def per_phase(self):
        """injection-free per-phase consumption (W, var) as seen by OpenDSS"""
        assert not self.delta()
        if self.kind == "sym_load":
            return [self.p / 3.0] * 3, [self.q / 3.0] * 3
        sign = -1.0 if self.kind == "asym_gen" else 1.0
        return [sign * x for x in self.p], [sign * x for x in self.q]

    def per_branch(self):
        """injection-free consumption (W, var) of the delta entries ab, bc and ca as seen by OpenDSS"""
        assert self.delta()
        sign = -1.0 if self.kind == "asym_gen" else 1.0
        return [sign * x for x in self.p], [sign * x for x in self.q]


def short_circuit_impedance(uk, pk, sn):
    """series impedance in per unit on sn from the short circuit voltage and the load losses, as PGM computes it"""
    r = pk / sn
    x_squared = uk * uk - r * r
    return complex(r, math.copysign(math.sqrt(x_squared) if x_squared > 0.0 else 0.0, uk))


WINDING_CONNECTION = {"a": 1, "b": 2, "c": 3, "ab": 4, "ac": 5, "bc": 6, "ab_center_tapped": 7}
# the centre-tapped (split-phase) to winding
CENTER_TAPPED = "ab_center_tapped"
BRANCH_SIDE = {"from": 0, "to": 1}


@dataclass
class Transformer1P:
    id: int
    from_node: int
    to_node: int
    # winding phases: "a", "b", "c" (phase-to-earth), "ab", "ac", "bc" (phase-to-phase) or, on the to side only,
    # "ab_center_tapped" (two halves, a to the centre tap and the centre tap to b)
    from_conn: str
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
    # centre-tapped to winding: the short circuit test between the primary and one half (on sn / 2); None means the
    # default (Kersting) split of uk and pk
    uk_half: float | None = None
    pk_half: float | None = None

    def center_tapped(self):
        return self.to_conn == CENTER_TAPPED

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
            "uk_half": self.uk_half,
            "pk_half": self.pk_half,
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

    def ct_impedances(self):
        """(z_f, z_hs, z_0, z_b) of a centre-tapped transformer in per unit on sn (the cases have no uk_min/max, so all
        at the nominal tap): z_f is the full winding (primary against both halves in series), z_hs the primary against
        one half (from uk_half/pk_half on sn / 2, or the Kersting default 1.5 r + j1.2 x), and z_0 = z_f - z_b / 2 and
        z_b = 2 (z_hs - z_f) the star branches of the primary and of each half"""
        z_f = short_circuit_impedance(self.uk, self.pk, self.sn)
        uk_half = 0.5 * abs(complex(1.5 * z_f.real, 1.2 * z_f.imag)) if self.uk_half is None else self.uk_half
        pk_half = 0.375 * self.pk if self.pk_half is None else self.pk_half
        z_hs = short_circuit_impedance(2.0 * uk_half, 4.0 * pk_half, self.sn)
        z_b = 2.0 * (z_hs - z_f)
        z_0 = z_f - 0.5 * z_b
        return z_f, z_hs, z_0, z_b

    def ct_primitive_admittance(self):
        """3x3 primitive admittance (siemens) of the primary and the two halves of a centre-tapped transformer, in the
        N-winding form A^T Z_B^-1 A with Z_B = z_0 1 + z_b I (per unit on sn), scaled with sn / (V_i V_j) on the tapped
        winding voltages u1', u2' / 2, u2' / 2; written independently of PGM, as the reference of the YPrim check"""
        _, _, z_0, z_b = self.ct_impedances()
        t1, t2 = self.taps()
        a = np.array([[1.0, -1.0, 0.0], [1.0, 0.0, -1.0]])
        z_branch = z_0 * np.ones((2, 2)) + z_b * np.eye(2)
        y_pu = a.T @ np.linalg.inv(z_branch) @ a
        v = np.array([self.u1 * t1, 0.5 * self.u2 * t2, 0.5 * self.u2 * t2])
        return y_pu * self.sn / np.outer(v, v)

    def dss_parameters(self):
        """the OpenDSS transformer definition without the buses"""
        if self.center_tapped():
            return self.dss_parameters_center_tapped()
        r_pct = 100.0 * self.pk / self.sn
        x_pct = 100.0 * math.sqrt(self.uk**2 - (self.pk / self.sn) ** 2)
        imag_pct = 100.0 * math.sqrt(self.i0**2 - (self.p0 / self.sn) ** 2)
        t1, t2 = self.taps()
        return (
            f"conns=[wye wye] kVs=[{self.u1 / 1e3!r} {self.u2 / 1e3!r}] kVAs=[{self.sn / 1e3!r} {self.sn / 1e3!r}]"
            f" XHL={x_pct!r} %loadloss={r_pct!r} %noloadloss={100.0 * self.p0 / self.sn!r} %imag={imag_pct!r}"
            f" ppm_antifloat=0 taps=[{t1!r} {t2!r}]"
        )

    def dss_parameters_center_tapped(self):
        """a 1-phase, 3-winding transformer (the primary and the two halves at u2 / 2), with equal kVAs = sn, the star
        resistances of z_0 and z_b per winding (%Rs, never %loadloss with three windings) and the reactances between
        the windings on sn: XHL = XHT = Im z_hs (primary against one half) and XLT = 2 Im z_b (half against half)"""
        _, z_hs, z_0, z_b = self.ct_impedances()
        assert self.i0 == 0.0, "the OpenDSS cases use i0 = p0 = 0, see the script"
        assert self.p0 == 0.0, "the OpenDSS cases use i0 = p0 = 0, see the script"
        t1, t2 = self.taps()
        kv_half = 0.5 * self.u2 / 1e3
        kva = self.sn / 1e3
        return (
            f"conns=[wye wye wye] kVs=[{self.u1 / 1e3!r} {kv_half!r} {kv_half!r}] kVAs=[{kva!r} {kva!r} {kva!r}]"
            f" %Rs=[{100.0 * z_0.real!r} {100.0 * z_b.real!r} {100.0 * z_b.real!r}]"
            f" XHL={100.0 * z_hs.imag!r} XHT={100.0 * z_hs.imag!r} XLT={200.0 * z_b.imag!r}"
            f" %noloadloss=0.0 %imag=0.0 ppm_antifloat=0 taps=[{t1!r} {t2!r} {t2!r}]"
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
            if conn == CENTER_TAPPED:
                assert side == "to", "a centre-tapped winding is only allowed on the to side"
                # half 1 from node 1 (slot a) to the centre tap, half 2 from the centre tap to node 2 (slot b); the
                # centre tap is earth (node 0) or node 4 with the electrode; clock 6 reverses both halves
                centre_tap = 0
                if r_g is not None or x_g is not None:
                    assert x_g is not None, "use x_grounding > 0 for an electrode"
                    assert x_g > 0.0, "use x_grounding > 0 for an electrode"
                    centre_tap = 4
                    shunts.append(
                        f"New Reactor.g{self.id}{side} phases=1 bus1={bus}.4 bus2={bus}.0 R={r_g or 0.0!r} X={x_g!r}"
                    )
                halves = [[1, centre_tap], [centre_tap, 2]]
                if self.clock == 6:
                    for half in halves:
                        half.reverse()
                buses += [bus + "".join(f".{t}" for t in half) for half in halves]
                continue
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
        transformer = (
            f"New Transformer.t{self.id} phases=1 windings={len(buses)} buses=[{' '.join(buses)}]"
            f" {self.dss_parameters()}"
        )
        return [transformer, *shunts]


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

    def has_bridge(self):
        """a bolted two-phase fault is a tiny resistance between two phases, whose accuracy is checked"""
        return self.fault_type == "two_phase" and self.r_f == 0.0

    def dss_elements(self, bridge_scale=1.0):
        """the OpenDSS elements of this fault, as (element, definition), e.g. ("Fault.f601", "phases=1 ..."); all
        element names start with f{id}; a fault with x_f is a Reactor (an OpenDSS Fault is resistive only); bridge_scale
        multiplies the resistance of a bolted two-phase fault (for the accuracy check)"""
        bus = f"n{self.node}"
        name = f"f{self.id}"
        if self.x_f != 0.0:
            return [(f"Reactor.{name}", self.dss_reactor(bus))]
        element = f"Fault.{name}"
        if self.fault_type == "single_phase_to_ground":
            k = "abc".index(self.fault_phase) + 1
            return [(element, f"phases=1 bus1={bus}.{k} r={self.dss_r()!r}")]
        if self.fault_type == "two_phase":
            k1, k2 = ("abc".index(p) + 1 for p in self.fault_phase)
            r = self.dss_r() * bridge_scale if self.has_bridge() else self.dss_r()
            return [(element, f"phases=1 bus1={bus}.{k1} bus2={bus}.{k2} r={r!r}")]
        if self.fault_type == "three_phase":
            # three phase-to-earth resistances r_f, as PGM models a three-phase fault with impedance
            return [(element, f"phases=3 bus1={bus}.1.2.3 bus2={bus}.0.0.0 r={self.dss_r()!r}")]
        if self.fault_type == "two_phase_to_ground":
            if self.r_f > 0.0:
                # PGM joins the two phases and connects them to earth through z_f; the bridge between the phases that
                # this needs in OpenDSS is too inaccurate, see the script docstring
                raise NotImplementedError("two_phase_to_ground with r_f > 0 has no accurate OpenDSS reference")
            # bolted: both phases to earth
            return [
                (f"{element}_{k}", f"phases=1 bus1={bus}.{k} r={DSS_BOLTED_FAULT_R!r}")
                for k in ("abc".index(p) + 1 for p in self.fault_phase)
            ]
        raise NotImplementedError(self.fault_type)

    def dss_reactor(self, bus):
        """the definition of the series R + jX Reactor of a fault with x_f, from the faulted phases to earth (or
        between the two phases of a two-phase fault); its R/X convention is checked by check_fault_impedance"""
        # OpenDSS may treat a Reactor with X = 0 specially, so the cases use x_f > 0
        assert self.x_f > 0.0, "use x_f > 0 for a fault with reactance"
        impedance = f"R={self.r_f!r} X={self.x_f!r}"
        if self.fault_type == "single_phase_to_ground":
            k = "abc".index(self.fault_phase) + 1
            return f"phases=1 bus1={bus}.{k} bus2={bus}.0 {impedance}"
        if self.fault_type == "two_phase":
            k1, k2 = ("abc".index(p) + 1 for p in self.fault_phase)
            return f"phases=1 bus1={bus}.{k1} bus2={bus}.{k2} {impedance}"
        if self.fault_type == "three_phase":
            return f"phases=3 bus1={bus}.1.2.3 bus2={bus}.0.0.0 {impedance}"
        raise NotImplementedError(f"{self.fault_type} with x_f has no OpenDSS reference")


@dataclass
class Case:
    path: str  # relative to tests/data
    u_rated: float
    nodes: list
    lines: list
    loads: list = field(default_factory=list)
    calculation: str = "power_flow"  # power_flow or short_circuit
    params: dict = field(default_factory=dict)
    # short circuit: per scenario a Fault or a tuple of simultaneous Faults; a batch if more than one scenario
    faults: list = field(default_factory=list)
    # cases that only check an error: (output type, raises name)
    raises_only: tuple | None = None
    source_id: int = 201
    sk: float = 100e6
    rx_ratio: float = 0.1
    z01_ratio: float = 1.0
    # u_rated per node where it differs from u_rated (the source node always has u_rated)
    node_u_rated: dict[int, float] = field(default_factory=dict)
    transformers: list = field(default_factory=list)
    # write the asym_load / asym_gen outputs (p, q, i per phase) to the reference as well
    load_output: bool = False

    def u_node(self, node):
        return self.node_u_rated.get(node, self.u_rated)


# ---------------------------------------------------------------------------------------------------------------
# PGM input


def fault_scenarios(case: Case) -> list[tuple[Fault, ...]]:
    """the faults of each scenario, as a tuple of simultaneous faults"""
    scenarios = [entry if isinstance(entry, tuple) else (entry,) for entry in case.faults]
    for scenario in scenarios:
        assert len({f.id for f in scenario}) == len(scenario), "a fault ID occurs once per scenario"
    return scenarios


def fault_ids(case: Case) -> list[int]:
    """the fault IDs of all scenarios, in the order of their first appearance"""
    return list(dict.fromkeys(f.id for scenario in fault_scenarios(case) for f in scenario))


def first_appearance(case: Case, fault_id: int) -> Fault:
    return next(f for scenario in fault_scenarios(case) for f in scenario if f.id == fault_id)


def voltage_factor(case: Case) -> float:
    """the IEC 60909 voltage factor c as PGM picks it, from the u_rated of the source node (not the fault location)"""
    if case.params["short_circuit_voltage_scaling"] == "maximum":
        return 1.1
    return 0.95 if case.u_node(case.nodes[0]) <= 1000.0 else 1.0


def pgm_input(case: Case, scenario: tuple[Fault, ...]):
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
    if case.faults:
        # every fault ID: the faults of the scenario with status 1, the others with status 0 (at their first appearance)
        present = {f.id: f for f in scenario}
        rows = []
        for fault_id in fault_ids(case):
            row = pgm_fault(present.get(fault_id, first_appearance(case, fault_id)), with_status=True)
            if fault_id not in present:
                row["status"] = 0
            rows.append(row)
        data["fault"] = rows
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


def fault_update(case: Case, scenario: tuple[Fault, ...]):
    """the fault update of a scenario: with a single fault ID the fault without status; otherwise every fault ID, the
    faults of the scenario with status 1 and the others switched off"""
    ids = fault_ids(case)
    if len(ids) == 1:
        return [pgm_fault(scenario[0], with_status=False)]
    present = {f.id: f for f in scenario}
    return [pgm_fault(present[i], with_status=True) if i in present else {"id": i, "status": 0} for i in ids]


# ---------------------------------------------------------------------------------------------------------------
# OpenDSS model


def set_active_element(name):
    """make an OpenDSS element the active one; OpenDSS keeps the previous element active for an unknown name, which
    would silently read the wrong element, so that is an error"""
    assert dss.Circuit.SetActiveElement(name) >= 0, f"OpenDSS has no element {name}"


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


def split_phase_nodes(case: Case):
    """the split-phase nodes, as PGM finds them (decision D25): the to nodes of the connected centre-tapped windings,
    spread along the connected lines without phase c (a node with phase c is never split-phase)"""
    split = {t.to_node for t in case.transformers if t.center_tapped() and t.from_status and t.to_status}
    changed = True
    while changed:
        changed = False
        for line in case.lines:
            if "c" in line.phases or not (line.from_status and line.to_status):
                continue
            for a, b in ((line.from_node, line.to_node), (line.to_node, line.from_node)):
                if a in split and b not in split:
                    split.add(b)
                    changed = True
    return split


def delta_kv(case: Case, node):
    """the OpenDSS kV of a delta load (the rated branch voltage, which PGM uses as U_delta,nom, decision D29):
    u_rated at a three-phase-pattern node and 2 u_rated / sqrt(3) (480 V on 240-0-240 V) at a split-phase node"""
    factor = 2.0 if node in split_phase_nodes(case) else math.sqrt(3.0)
    return factor * case.u_node(node) / math.sqrt(3.0) / 1e3


def dss_load_elements(case: Case, load: Load):
    """the OpenDSS Load elements of an active load, as (name, definition, branches); all names start with ld{id}.
    A wye load is 1-phase wye loads per phase with kV = u_rated / sqrt(3); a delta load is a 1-phase conn=delta load
    per non-zero entry (bus1 = n.p.q, kV = the branch voltage) or, with dss_three_phase, one native 3-phase conn=delta
    load with the total power. branches lists (p, q, S in VA) per delta entry of the element (for the base check)."""
    if not load.status:
        return []
    model = DSS_LOAD_MODEL[load.type]
    elements: list[tuple[str, str, list[tuple[int, int, complex]] | None]] = []
    if not load.delta():
        kv_phase = case.u_node(load.node) / math.sqrt(3.0) / 1e3
        p, q = load.per_phase()
        for k in range(3):
            if p[k] == 0.0 and q[k] == 0.0:
                continue
            definition = (
                f"bus1=n{load.node}.{k + 1} phases=1 conn=wye model={model} kV={kv_phase!r} kW={p[k] / 1e3!r}"
                f" kvar={q[k] / 1e3!r} vminpu=0.1 vmaxpu=2"
            )
            elements.append((f"ld{load.id}_{'abc'[k]}", definition, None))
        return elements
    kv = delta_kv(case, load.node)
    p, q = load.per_branch()
    if load.dss_three_phase:
        assert p[0] == p[1] == p[2], "a native 3-phase delta load has equal entries"
        assert q[0] == q[1] == q[2], "a native 3-phase delta load has equal entries"
        assert load.node not in split_phase_nodes(case)
        definition = (
            f"bus1=n{load.node}.1.2.3 phases=3 conn=delta model={model} kV={kv!r} kW={3.0 * p[0] / 1e3!r}"
            f" kvar={3.0 * q[0] / 1e3!r} vminpu=0.1 vmaxpu=2"
        )
        branches = [(a, b, complex(p[k], q[k])) for k, (a, b) in enumerate(DELTA_BRANCHES)]
        return [(f"ld{load.id}", definition, branches)]
    for k, (a, b) in enumerate(DELTA_BRANCHES):
        if p[k] == 0.0 and q[k] == 0.0:
            continue
        definition = (
            f"bus1=n{load.node}.{a}.{b} phases=1 conn=delta model={model} kV={kv!r} kW={p[k] / 1e3!r}"
            f" kvar={q[k] / 1e3!r} vminpu=0.1 vmaxpu=2"
        )
        elements.append((f"ld{load.id}_{DELTA_NAMES[k]}", definition, [(a, b, complex(p[k], q[k]))]))
    return elements


def dss_script(case: Case, faults: tuple[Fault, ...], voltage_pu: float, with_loads: bool, bridge_scale=1.0):
    """the OpenDSS script of a power flow (no faults) or of a short-circuit scenario (its simultaneous faults)"""
    z1, z0 = source_z(case)
    if faults:
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
            lines += [f"New Load.{name} {definition}" for name, definition, _ in dss_load_elements(case, load)]
    for fault in faults:
        lines += [f"New {element} {definition}" for element, definition in fault.dss_elements(bridge_scale)]
    bases = ", ".join(repr(u / 1e3) for u in sorted({case.u_node(n) for n in case.nodes}, reverse=True))
    lines += [f"Set VoltageBases=[{bases}]", "CalcVoltageBases", "Set Tolerance=1e-12", "Set MaxIterations=200"]
    if faults:
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
        set_active_element(f"Line.l{line.id}")
        assert dss.CktElement.NumPhases() == n, f"OpenDSS line {line.id} has {dss.CktElement.NumPhases()} phases"
        assert dss.CktElement.NumConductors() == n, f"OpenDSS line {line.id} still has a neutral conductor"
        assert dss.CktElement.NumTerminals() == 2
        assert dss.Lines.Phases() == n
        dss.LineCodes.Name(f"lc{line.id}")
        assert dss.LineCodes.Phases() == n, f"OpenDSS did not Kron-reduce the LineCode of line {line.id}"
        set_active_element(f"Line.l{line.id}")
        r = np.array(dss.Lines.RMatrix()).reshape(n, n) * line.length
        x = np.array(dss.Lines.XMatrix()).reshape(n, n) * line.length
        expected = line.z_phase()
        dev = np.max(np.abs(r + 1j * x - expected))
        assert dev < 1e-12 * np.max(np.abs(expected)), f"OpenDSS Kron differs from numpy for line {line.id}: {dev}"
        print(f"  OpenDSS-reduced line {line.id}: {n} phases, matches numpy Kron (max dev {dev:.1e} ohm)")


def bus_phasors(node):
    """complex voltages (V) of the phase nodes 1 to 3 of a bus, by node number"""
    dss.Circuit.SetActiveBus(f"n{node}")
    raw = np.array(dss.Bus.Voltages()).reshape(-1, 2)
    return {number: complex(*pair) for number, pair in zip(dss.Bus.Nodes(), raw) if 1 <= number <= 3}


def check_delta_load_bases(case: Case):
    """the decisive check of the OpenDSS delta-load conventions, after the solve: the total power of every delta Load
    element (the sum of CktElement.Powers over its conductors) equals sum_k S_k (|V_p - V_q| / kV)^n with the OpenDSS
    voltages, S_k the entry, kV the branch voltage of the element and n = 0, 1, 2 for model 1, 5, 2, to 1e-9 relative.
    This proves that the kV of a delta load (1-phase and native 3-phase) is the branch voltage and that model 5 and 2
    scale with the branch voltage as PGM does (D29). It compares the element total only, which for const_pq equals 3S
    whatever the split, so it does not prove the split of a native 3-phase delta load into a third per branch; that is
    proven end-to-end by the per-phase load outputs matching the OpenDSS terminal values"""
    exponent = {"const_pq": 0, "const_i": 1, "const_y": 2}
    worst = 0.0
    n_elements = 0
    for load in case.loads:
        if not load.delta():
            continue
        kv = delta_kv(case, load.node)
        v = bus_phasors(load.node)
        for name, _, branches in dss_load_elements(case, load):
            set_active_element(f"Load.{name}")
            raw = np.array(dss.CktElement.Powers()).reshape(-1, 2)
            total = complex(*raw.sum(axis=0)) * 1e3
            expected = sum(s * (abs(v[a] - v[b]) / (kv * 1e3)) ** exponent[load.type] for a, b, s in branches)
            dev = abs(total - expected) / abs(expected)
            assert dev < 1e-9, f"OpenDSS delta load {name}: {total} VA instead of {expected} VA ({dev:.1e} relative)"
            worst = max(worst, dev)
            n_elements += 1
    if n_elements:
        print(
            f"  delta load bases: {n_elements} OpenDSS delta loads agree with the branch-voltage model"
            f" (max dev {worst:.1e} relative)"
        )


def check_fault_bridge(case: Case, scenario: tuple[Fault, ...], currents):
    """a bolted two-phase fault (a bridge of 1e-13 ohm between the phases) solved again with a 10x larger bridge gives
    the same fault current to 1e-10 relative: the reference is limited neither by the bridge resistance nor by the
    round-off of its large conductance; currents are the fault currents of the scenario, per fault"""
    if not any(fault.has_bridge() for fault in scenario):
        return
    run_dss(dss_script(case, scenario, voltage_factor(case), with_loads=False, bridge_scale=10.0))
    for fault, i_f in zip(scenario, currents):
        if not fault.has_bridge():
            continue
        dev = np.max(np.abs(fault_currents(fault) - i_f)) / np.max(np.abs(i_f))
        assert dev < 1e-10, (
            f"bolted two-phase fault at node {fault.node} is not accurate in OpenDSS: {dev:.1e} relative"
        )
        print(f"  bolted two-phase fault at node {fault.node}: a 10x larger bridge changes i_f by {dev:.1e} relative")


def check_fault_impedance(fault: Fault, i_f):
    """the decisive check of the R/X convention of the Reactor of a fault with x_f, after the solve: with the OpenDSS
    node voltages, its current equals V / (r_f + j x_f) to 1e-9 relative on each faulted phase, with V the phase
    voltage to earth, or V_k1 - V_k2 (and -(V_k1 - V_k2) on k2) for a two-phase fault"""
    if fault.x_f == 0.0:
        return
    z_f = complex(fault.r_f, fault.x_f)
    v = bus_phasors(fault.node)
    expected = np.zeros(3, dtype=complex)
    if fault.fault_type == "two_phase":
        k1, k2 = ("abc".index(p) + 1 for p in fault.fault_phase)
        expected[k1 - 1] = (v[k1] - v[k2]) / z_f
        expected[k2 - 1] = -expected[k1 - 1]
    else:
        for p in fault.fault_phase:
            k = "abc".index(p) + 1
            expected[k - 1] = v[k] / z_f
    dev = np.max(np.abs(i_f - expected)) / np.max(np.abs(expected))
    assert dev < 1e-9, f"fault {fault.id} at node {fault.node}: the Reactor current is not V / z_f ({dev:.1e} relative)"
    print(
        f"  fault {fault.id} at node {fault.node}: the Reactor current equals V / z_f with z_f = {z_f} ohm"
        f" (dev {dev:.1e} relative)"
    )


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
        set_active_element("Transformer.t")
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


def check_center_tapped_transformer(transformer: Transformer1P):
    """sanity checks of the OpenDSS centre-tapped (3-winding) transformer on its own, fed with the from-side winding
    voltage:

    - YPrim: the 6x6 primitive admittance of OpenDSS (siemens; conductors [p+, p-, h1+, h1-, h2+, h2-]) equals
      C^T Y C of the numpy N-winding model (C has the rows e_p+ - e_p-, e_h1+ - e_h1-, e_h2+ - e_h2-); this is the
      decisive check of the three-winding conventions (%Rs, XHL/XHT/XLT on kVAs[1], kVs and taps per winding)
    - open circuit: U_b = -U_a and |U_a| = u2 / 2 (times the tap ratio)
    """
    assert transformer.tap_side in (None, "from"), "a to-side tap changes the ohmic impedance differently"
    t1, t2 = transformer.taps()
    u1 = transformer.u1
    dss.Basic.ClearAll()
    for command in (
        "Clear",
        f"Set DefaultBaseFrequency={FREQ:g}",
        f"New Circuit.t phases=1 basekv={u1 / 1e3!r} pu=1 angle=0 bus1=src.1 Z1=[0, 1e-9] Z0=[0, 1e-9]",
        f"New Transformer.t phases=1 windings=3 buses=[src.1.0 sec.1.0 sec.0.2] {transformer.dss_parameters()}",
        "Set VoltageBases=[1]",
        "CalcVoltageBases",
        "Set Tolerance=1e-12",
        "Solve",
    ):
        dss.Text.Command(command)
    assert dss.Solution.Converged()
    set_active_element("Transformer.t")
    assert dss.CktElement.NumTerminals() == 3
    assert dss.CktElement.NumConductors() == 2
    raw = np.array(dss.CktElement.YPrim())
    y_dss = (raw[0::2] + 1j * raw[1::2]).reshape(6, 6)
    c = np.zeros((3, 6))
    for w in range(3):
        c[w, 2 * w] = 1.0
        c[w, 2 * w + 1] = -1.0
    y_expected = c.T @ transformer.ct_primitive_admittance() @ c
    dev = np.max(np.abs(y_dss - y_expected)) / np.max(np.abs(y_expected))
    assert dev < 1e-9, f"YPrim of centre-tapped transformer {transformer.id} differs by {dev:.1e} relative"

    dss.Circuit.SetActiveBus("src")
    v1 = complex(*dss.Bus.Voltages()[:2])
    assert abs(abs(v1) / u1 - 1.0) < 1e-9, f"the 1-phase source is not at the winding voltage: {abs(v1)}"
    dss.Circuit.SetActiveBus("sec")
    v = dict(zip(dss.Bus.Nodes(), (complex(*pair) for pair in np.array(dss.Bus.Voltages()).reshape(-1, 2))))
    ratio = v[1] / v1
    expected_ratio = 0.5 * transformer.u2 * t2 / (u1 * t1)
    dev_ratio = abs(ratio / expected_ratio - 1.0)
    dev_opposite = abs(v[2] + v[1]) / abs(v[1])
    assert dev_ratio < 1e-9, f"open-circuit ratio {ratio} instead of {expected_ratio}"
    assert dev_opposite < 1e-9, f"open circuit: U_b = {v[2]} is not -U_a = {-v[1]}"
    print(
        f"  centre-tapped transformer {transformer.id}: YPrim agrees (max dev {dev:.1e} relative); open circuit"
        f" |U_a| = u2 / 2 (dev {dev_ratio:.1e}) and U_b = -U_a (dev {dev_opposite:.1e})"
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
    """per side, per phase complex values of the active element: the first terminal is the from side, and the other
    terminals (one for a line or a two-winding transformer; the two halves of a centre-tapped transformer) are summed
    per node into the to side; nodes 0 and 4 are dropped"""
    set_active_element(name)
    raw = values_getter()
    order = dss.CktElement.NodeOrder()
    n_terminals = dss.CktElement.NumTerminals()
    result = [np.zeros(3, dtype=complex), np.zeros(3, dtype=complex)]
    assigned = [set(), set()]
    for terminal in range(n_terminals):
        side = min(terminal, 1)
        for k in range(n_cond):
            idx = terminal * n_cond + k
            number = order[idx]
            if 1 <= number <= 3:
                value = complex(raw[2 * idx], raw[2 * idx + 1])
                # assign the first value (keeping the sign of a zero), add the next ones
                if number in assigned[side]:
                    result[side][number - 1] += value
                else:
                    result[side][number - 1] = value
                    assigned[side].add(number)
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
    if case.load_output:
        assert not isolated, "the load reference assumes energised loads"
        for kind in ("asym_load", "asym_gen"):
            rows = [load_pf_reference(case, load) for load in case.loads if load.kind == kind]
            if rows:
                result[kind] = rows
        assert all(load.kind != "sym_load" for load in case.loads), "no sym_load reference"
    return result


def load_pf_reference(case: Case, load: Load):
    """per-phase terminal values of a load or generator: the sums per node over its OpenDSS Load elements of the
    powers (p, q; negated for a generator) and of the currents (i = the magnitude of the sum); for a delta load these
    are the phase quantities that PGM reports (decision D30)"""
    s = np.zeros(3, dtype=complex)
    current = np.zeros(3, dtype=complex)
    for name, _, _ in dss_load_elements(case, load):
        set_active_element(f"Load.{name}")
        n_cond = dss.CktElement.NumConductors()
        s += terminal_values(f"Load.{name}", dss.CktElement.Powers, n_cond)[0]
        current += terminal_values(f"Load.{name}", dss.CktElement.Currents, n_cond)[0]
    sign = -1.0 if load.kind == "asym_gen" else 1.0
    s = sign * s * 1e3
    return {
        "id": load.id,
        "p": [float(v.real) for v in s],
        "q": [float(v.imag) for v in s],
        "i": [float(abs(v)) for v in current],
    }


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
    # each half of a centre-tapped winding is rated sn / 2
    to_power = 2.0 * max(abs(s_to[0]), abs(s_to[1])) if transformer.center_tapped() else abs(np.sum(s_to))
    return {
        "id": transformer.id,
        "p_from": [float(v.real) for v in s_from],
        "q_from": [float(v.imag) for v in s_from],
        "p_to": [float(v.real) for v in s_to],
        "q_to": [float(v.imag) for v in s_to],
        "loading": float(max(abs(np.sum(s_from)), to_power) / transformer.sn),
    }


def sc_reference(case: Case, scenario: tuple[Fault, ...], currents):
    """the reference of a short-circuit scenario; currents are the fault currents of the scenario, per fault; a fault
    ID that is not in the scenario (status 0) gets i_f = 0, as PGM's null output"""
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
    result = {"node": nodes, "asym_line": lines}
    if case.transformers:
        # the fault current is only printed with the electrode current, for a single fault
        i_f = currents[0] if len(scenario) == 1 else None
        result["single_phase_transformer"] = [transformer_sc_reference(t, i_f) for t in case.transformers]
    present = {fault.id: i_f for fault, i_f in zip(scenario, currents)}
    rows = []
    for fault_id in fault_ids(case):
        i_f_mag, i_f_angle = mag_angle(present.get(fault_id, np.zeros(3, dtype=complex)))
        rows.append({"id": fault_id, "i_f": i_f_mag, "i_f_angle": i_f_angle})
    result["fault"] = rows
    return result


def fault_currents(fault: Fault):
    """the current from the node into the fault on each faulted phase: the sum over all OpenDSS elements of the fault
    and over all their terminals (a terminal at earth, node 0, is dropped)"""
    i_f = np.zeros(3, dtype=complex)
    assigned = set()
    for element, _ in fault.dss_elements():
        set_active_element(element)
        raw = dss.CktElement.Currents()
        order = dss.CktElement.NodeOrder()
        for idx, number in enumerate(order):
            if 1 <= number <= 3:
                value = complex(raw[2 * idx], raw[2 * idx + 1])
                # assign the first value (keeping the sign of a zero), add the next ones
                if number in assigned:
                    i_f[number - 1] += value
                else:
                    i_f[number - 1] = value
                    assigned.add(number)
    return i_f


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


def check_electrodes(transformer: Transformer1P, i_f=None):
    """a consistency check of the OpenDSS netlist (KCL), not of PGM: the current from node 4 into the electrode equals
    the current into the winding at its phase terminal (the winding has no other path), so the winding current that
    PGM is compared with is also the electrode current; for a centre-tapped winding, it equals the sum I_a + I_b of
    the currents into the two halves at the legs (KCL at node 4, the centre tap); i_f is the fault current, if any"""
    for side, terminal, r_g, x_g, status in (
        ("from", 0, transformer.r_g_from, transformer.x_g_from, transformer.from_status),
        ("to", 1, transformer.r_g_to, transformer.x_g_to, transformer.to_status),
    ):
        if (r_g is None and x_g is None) or not status:
            continue
        if side == "to" and transformer.center_tapped():
            check_center_tap_electrode(transformer, i_f)
            continue
        set_active_element(f"Transformer.t{transformer.id}")
        raw = dss.CktElement.Currents()
        order = dss.CktElement.NodeOrder()
        winding = {order[idx]: complex(raw[2 * idx], raw[2 * idx + 1]) for idx in range(2 * terminal, 2 * terminal + 2)}
        phase = "abc".index(transformer.from_conn if side == "from" else transformer.to_conn) + 1
        set_active_element(f"Reactor.g{transformer.id}{side}")
        i_electrode = complex(*dss.CktElement.Currents()[:2])
        dev = abs(i_electrode - winding[phase])
        assert dev < 1e-9 * abs(i_electrode) + 1e-9, f"electrode of transformer {transformer.id} {side}: {dev}"
        assert abs(winding[4] + winding[phase]) < 1e-9 * abs(i_electrode) + 1e-9
        fault_text = "" if i_f is None else f"; fault current on phase {'abc'[phase - 1]}: {abs(i_f[phase - 1]):.6g} A"
        print(
            f"  electrode of transformer {transformer.id} ({side}): {abs(i_electrode):.6g} A, equal to the winding"
            f" current (dev {dev:.1e} A){fault_text}"
        )


def check_center_tap_electrode(transformer: Transformer1P, i_f=None):
    """KCL at node 4, the centre tap: the electrode current equals I_a + I_b, the currents into the two halves at the
    legs (nodes 1 and 2; the terminals 2 and 3 of the transformer); i_f is the fault current, if any"""
    set_active_element(f"Transformer.t{transformer.id}")
    raw = dss.CktElement.Currents()
    order = dss.CktElement.NodeOrder()
    legs = 0j
    centre = 0j
    for idx in range(2, 6):
        value = complex(raw[2 * idx], raw[2 * idx + 1])
        if order[idx] in (1, 2):
            legs += value
        else:
            assert order[idx] == 4
            centre += value
    set_active_element(f"Reactor.g{transformer.id}to")
    i_electrode = complex(*dss.CktElement.Currents()[:2])
    dev = abs(i_electrode - legs)
    assert dev < 1e-9 * abs(i_electrode) + 1e-9, f"centre-tap electrode of transformer {transformer.id}: {dev}"
    assert abs(centre + legs) < 1e-9 * abs(i_electrode) + 1e-9
    fault_text = "" if i_f is None else f"; fault current per phase: {', '.join(f'{abs(v):.6g}' for v in i_f)} A"
    print(
        f"  centre-tap electrode of transformer {transformer.id}: {abs(i_electrode):.6g} A, equal to I_a + I_b"
        f" (dev {dev:.1e} A){fault_text}"
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
    scenarios = fault_scenarios(case)
    write_json(target / "input.json", pgm_input(case, scenarios[0] if scenarios else ()))
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
        if transformer.center_tapped():
            check_center_tapped_transformer(transformer)
        else:
            check_transformer(transformer)

    if case.calculation == "power_flow":
        script = dss_script(case, (), 1.0, with_loads=True)
        run_dss(script, case)
        check_delta_load_bases(case)
        for transformer in case.transformers:
            if transformer.center_tapped():
                check_electrodes(transformer)
        write_text(target / "test.dss", "\n".join(script) + "\n")
        write_json(target / "asym_output.json", output_dataset("asym_output", pf_reference(case)))
        return

    scripts = []
    references = []
    for scenario in scenarios:
        script = dss_script(case, scenario, voltage_factor(case), with_loads=False)
        run_dss(script, case)
        scripts.append(script)
        currents = [fault_currents(fault) for fault in scenario]
        for fault, i_f in zip(scenario, currents):
            check_fault_impedance(fault, i_f)
        references.append(sc_reference(case, scenario, currents))
        check_fault_bridge(case, scenario, currents)
    dss_text = []
    for k, script in enumerate(scripts):
        if len(scripts) > 1:
            dss_text.append(f"! scenario {k}")
        dss_text += script
        dss_text.append("")
    write_text(target / "test.dss", "\n".join(dss_text))
    if len(scenarios) == 1:
        write_json(target / "sc_output.json", output_dataset("sc_output", references[0]))
    else:
        write_json(
            target / "update_batch.json",
            {
                "version": "1.0",
                "type": "update",
                "is_batch": True,
                "attributes": {},
                "data": [{"fault": fault_update(case, scenario)} for scenario in scenarios],
            },
        )
        write_json(target / "sc_output_batch.json", output_dataset("sc_output", references, is_batch=True))


# ---------------------------------------------------------------------------------------------------------------
# cases

# Tolerances: per case, the smallest round rtol (with atol <= 1e-6), then the smallest round atol, for which the worst
# |actual - expected| / (atol + rtol |expected|) is below 0.1 (a 10x margin); see reduced_phase_compare_pgm.py.
# PGM and OpenDSS agree to about 3.4e-10 relative everywhere. iterative_current runs with IC_ERROR_TOLERANCE (below), so
# it is no longer limited by its convergence criterion and every method agrees about equally well; the cases that need
# more than (1e-8, 1e-8) need it for newton_raphson as well: near-zero line powers (reduced-phase-feeder, swer-a and
# swer-chain: atol) and the nearly cancelling reactive power of the lines (delta two-phase: rtol).
SC_PARAMS = {
    "calculation_method": "iec60909",
    "rtol": 1e-8,
    "atol": 1e-8,
    "short_circuit_voltage_scaling": "maximum",
    "experimental_features": "enabled",
}


def sc_raises_params(raises, reason):
    return dict(SC_PARAMS, raises={"raises": raises, "reason": reason})


# error_tolerance of iterative_current only (extra_params), two orders below the default 1e-8: at 1e-8 iterative_current
# stops on its convergence residual (about error_tolerance * rho / (1 - rho) on the voltage), e.g. both
# split-phase-flat-start cases missed (1e-8, 1e-8) with a worst ratio of up to 2.645; at 1e-10 it agrees with OpenDSS as
# closely as newton_raphson does (ratio 0.034 there, in 8 and 9 iterations; at most 11 iterations in any case,
# phase-to-phase). It is per method, because newton_raphson stalls at about 5e-11 on the M10 networks and does not
# converge for an error_tolerance of 1e-11 or less.
IC_ERROR_TOLERANCE = 1e-10


def pf_params(methods, rtol, atol):
    params = {"calculation_method": methods, "rtol": rtol, "atol": atol, "experimental_features": "enabled"}
    if "iterative_current" in methods:
        params["extra_params"] = {"iterative_current": {"error_tolerance": IC_ERROR_TOLERANCE}}
    return params


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
        params=pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-6),
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


def swer_short_circuit_cases():
    """the short-circuit phase validity cases (M8) in the swer-chain network without loads: 22 kV abc nodes 1-2, the
    isolating transformer ab -> a with its electrode, SWER nodes 3-4 and 240 V nodes 5 (behind the tapped
    transformer) and 6"""
    base = swer_chain_short_circuit_case()
    spg = "single_phase_to_ground"
    return [
        # x_f != 0: an OpenDSS Reactor, checked by check_fault_impedance
        replace(
            base,
            path="short_circuit/swer-fault-impedance",
            faults=[
                # the fault current returns through the electrode of the isolating transformer
                Fault(601, 4, spg, "a", 10.0, 5.0),
                Fault(601, 6, spg, "a", 0.05, 0.02),
                Fault(601, 4, spg, "a", 0.0, 8.0),
            ],
        ),
        # simultaneous faults; every fault has an impedance that keeps the voltage of the faults downstream of it, so
        # that every fault current is far from zero
        replace(
            base,
            path="short_circuit/swer-multiple-faults",
            faults=[
                (Fault(601, 4, spg, "a", 10.0), Fault(602, 6, spg, "a", 0.05)),
                # a three-phase node and SWER and 240 V nodes in one calculation
                (Fault(601, 2, spg, "a", 5.0), Fault(602, 5, spg, "a", 0.0), Fault(603, 3, spg, "a", 20.0, 10.0)),
            ],
        ),
        # c = 1.0 from the 22 kV source node, also for the 240 V fault (IEC 60909: 0.95 at the fault location)
        replace(
            base,
            path="short_circuit/swer-c-minimum",
            params=dict(SC_PARAMS, short_circuit_voltage_scaling="minimum"),
            faults=[Fault(601, 4, spg, "a", 0.0), Fault(601, 5, spg, "a", 0.0)],
        ),
        replace(
            base,
            path="short_circuit/swer-fault-on-absent-phase",
            faults=[Fault(601, 4, spg, "b", 0.0)],
            params=sc_raises_params(
                "InvalidShortCircuitPhaseOrType",
                "The short circuit phases (b) of fault 601 are not all present at its node 4 (node phases: a); this"
                " short circuit type cannot be applied there",
            ),
            raises_only=("sc_output", "InvalidShortCircuitPhaseOrType"),
        ),
    ]


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


def split_phase_case(path, all_const_y, params):
    """22 kV abc -> SWER isolating transformer ab -> a (electrode) -> SWER line -> a -> 240-0-240 V centre-tapped
    (Kersting defaults, solid centre tap, from-side tap) with an ab service line; and 22 kV bc -> 240-0-240 V
    centre-tapped (clock 6, uk_half/pk_half given, centre-tap electrode)"""

    def t(load_type):
        return "const_y" if all_const_y else load_type

    return Case(
        path=path,
        u_rated=SWER_U,
        nodes=[1, 2, 3, 4, 5, 6, 7],
        node_u_rated={5: LV_U, 6: LV_U, 7: LV_U},
        lines=[
            Line(101, 1, 2, "abc", 2.0),
            Line(102, 3, 4, "a", 20.0, swer=True),
            # the service line off the split-phase node, which carries the split-phase slot pattern
            Line(103, 5, 7, "ab", 0.3, neutral=True),
        ],
        transformers=[
            Transformer1P(151, 2, 3, "ab", "a", 22e3, 12.7e3, 200e3, 0.05, 2000.0, r_g_to=5.0, x_g_to=0.1),
            Transformer1P(
                152,
                4,
                5,
                "a",
                CENTER_TAPPED,
                12.7e3,
                480.0,
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
            Transformer1P(
                153,
                2,
                6,
                "bc",
                CENTER_TAPPED,
                22e3,
                480.0,
                50e3,
                0.035,
                500.0,
                clock=6,
                r_g_to=0.5,
                x_g_to=0.05,
                uk_half=0.022,
                pk_half=200.0,
            ),
        ],
        # wye loads on the legs a and b (240 V to true earth), unbalanced
        loads=[
            Load(301, 5, "asym_load", t("const_pq"), [8e3, 3e3, 0.0], [2.5e3, 1e3, 0.0]),
            Load(302, 5, "asym_load", "const_y", [1.5e3, 2.5e3, 0.0], [0.5e3, 0.8e3, 0.0]),
            Load(303, 7, "asym_load", t("const_pq"), [3e3, 1.5e3, 0.0], [1e3, 0.4e3, 0.0]),
            Load(304, 7, "asym_load", t("const_i"), [1e3, 2e3, 0.0], [0.3e3, 0.6e3, 0.0]),
            Load(305, 6, "asym_load", "const_y", [6e3, 2e3, 0.0], [2e3, 0.5e3, 0.0]),
            Load(306, 6, "asym_load", t("const_pq"), [1e3, 0.5e3, 0.0], [0.5e3, 0.2e3, 0.0]),
        ],
        params=params,
    )


def split_phase_short_circuit_case():
    """the split-phase network without loads; faults on the legs of both split-phase nodes and of the service line"""
    return replace(
        split_phase_case("short_circuit/single-phase-transformer-split-phase", all_const_y=False, params=SC_PARAMS),
        loads=[],
        calculation="short_circuit",
        faults=[
            # leg to earth (the solid centre tap) through half 1 or half 2
            Fault(601, 5, "single_phase_to_ground", "a", 0.0),
            Fault(601, 5, "single_phase_to_ground", "b", 0.1),
            # leg to leg (480 V), through both halves in series
            Fault(601, 5, "two_phase", "ab", 0.0),
            Fault(601, 5, "two_phase_to_ground", "ab", 0.0),
            # not degenerate (at the split-phase node it equals the two-phase fault): two phases of the 22 kV node
            Fault(601, 2, "two_phase_to_ground", "bc", 0.0),
            # at the end of the ab service line
            Fault(601, 7, "single_phase_to_ground", "a", 0.0),
            # the fault current returns through the centre-tap electrode of transformer 153
            Fault(601, 6, "single_phase_to_ground", "a", 0.0),
            Fault(601, 6, "two_phase", "ab", 0.05),
        ],
    )


DELTA_PATH = "power_flow/asym-load-delta/"


def delta_three_phase_case(path, all_const_y, params):
    """11 kV, abc only: delta loads of every type (a native 3-phase OpenDSS delta load and 1-phase delta loads), a
    wye load and a delta generator"""

    def t(load_type):
        return "const_y" if all_const_y else load_type

    return Case(
        path=path,
        u_rated=11e3,
        nodes=[1, 2, 3],
        lines=[Line(101, 1, 2, "abc", 2.0), Line(102, 2, 3, "abc", 1.5, neutral=True, c_form="c0c1")],
        loads=[
            Load(301, 2, "asym_load", t("const_pq"), [150e3] * 3, [50e3] * 3, connection="delta", dss_three_phase=True),
            Load(302, 3, "asym_load", t("const_i"), [120e3, 0.0, 80e3], [40e3, 0.0, 20e3], connection="delta"),
            Load(303, 3, "asym_load", "const_y", [0.0, 100e3, 60e3], [0.0, 30e3, 15e3], connection="delta"),
            Load(304, 3, "asym_load", t("const_pq"), [60e3, 30e3, 0.0], [20e3, 10e3, 0.0]),
            Load(401, 2, "asym_gen", t("const_pq"), [0.0, 50e3, 0.0], [0.0, 10e3, 0.0], connection="delta"),
        ],
        params=params,
        load_output=True,
    )


def delta_two_phase_case(path, all_const_y, params):
    """11 kV abc -> ab laterals (with and without neutral) and an ac lateral: phase-to-phase loads on two-phase nodes
    (ab and ca entries), a wye load on a two-phase node and a delta load on the abc node"""

    def t(load_type):
        return "const_y" if all_const_y else load_type

    return Case(
        path=path,
        u_rated=11e3,
        nodes=[1, 2, 3, 4, 5],
        lines=[
            Line(101, 1, 2, "abc", 2.0),
            Line(102, 2, 3, "ab", 2.0, neutral=True),
            Line(103, 3, 4, "ab", 1.0),
            Line(104, 2, 5, "ac", 1.5),
        ],
        loads=[
            Load(301, 3, "asym_load", t("const_pq"), [80e3, 0.0, 0.0], [25e3, 0.0, 0.0], connection="delta"),
            Load(302, 4, "asym_load", t("const_i"), [50e3, 0.0, 0.0], [15e3, 0.0, 0.0], connection="delta"),
            Load(303, 4, "asym_load", "const_y", [20e3, 10e3, 0.0], [5e3, 3e3, 0.0]),
            Load(304, 5, "asym_load", t("const_pq"), [0.0, 0.0, 60e3], [0.0, 0.0, 20e3], connection="delta"),
            Load(305, 2, "asym_load", "const_y", [40e3, 30e3, 20e3], [10e3, 8e3, 5e3], connection="delta"),
        ],
        params=params,
        load_output=True,
    )


def delta_split_phase_case(path, all_const_y, params):
    """the split-phase network of M6 with 480 V (leg-to-leg) delta loads added at both split-phase nodes and at the
    end of the service line"""

    def t(load_type):
        return "const_y" if all_const_y else load_type

    case = split_phase_case(path, all_const_y, params)
    return replace(
        case,
        loads=[
            *case.loads,
            Load(307, 5, "asym_load", t("const_pq"), [4e3, 0.0, 0.0], [1e3, 0.0, 0.0], connection="delta"),
            Load(308, 7, "asym_load", t("const_i"), [2e3, 0.0, 0.0], [0.5e3, 0.0, 0.0], connection="delta"),
            Load(309, 6, "asym_load", "const_y", [5e3, 0.0, 0.0], [1e3, 0.0, 0.0], connection="delta"),
        ],
        load_output=True,
    )


def raises_params(methods, raises, reason):
    return {
        "calculation_method": methods,
        "rtol": 1e-6,
        "atol": 1e-4,
        "experimental_features": "enabled",
        "raises": {"raises": raises, "reason": reason},
    }


def small_delta_case(path, raises_only, params):
    """a small abc network with a delta load"""
    return Case(
        path=path,
        u_rated=10e3,
        nodes=[1, 2],
        lines=[Line(101, 1, 2, "abc", 5.0)],
        loads=[Load(301, 2, "asym_load", "const_pq", [200e3, 0.0, 0.0], [50e3, 0.0, 0.0], connection="delta")],
        params=params,
        raises_only=raises_only,
    )


def delta_cases():
    """the delta-load cases (M7); the tolerances follow the rule above: all methods agree to about 3.4e-10 relative
    (iterative_current with IC_ERROR_TOLERANCE), except that two-phase and two-phase-linear are limited by the
    reactive power of the lines, which nearly cancels (6e-6 var on a few kvar)"""
    return [
        delta_three_phase_case(
            DELTA_PATH + "three-phase",
            all_const_y=False,
            params=pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-8),
        ),
        delta_three_phase_case(
            DELTA_PATH + "three-phase-linear",
            all_const_y=True,
            params=pf_params(["linear", "linear_current"], 1e-8, 1e-8),
        ),
        delta_two_phase_case(
            DELTA_PATH + "two-phase",
            all_const_y=False,
            params=pf_params(["newton_raphson", "iterative_current"], 1e-7, 1e-8),
        ),
        delta_two_phase_case(
            DELTA_PATH + "two-phase-linear",
            all_const_y=True,
            params=pf_params(["linear", "linear_current"], 1e-7, 1e-8),
        ),
        delta_split_phase_case(
            DELTA_PATH + "split-phase",
            all_const_y=False,
            params=pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-8),
        ),
        # iterative_current ignores calculation_initialization and always uses the average-source start, which
        # equals the flat start here (one source, u_ref 1.0); the flat start itself is exercised by newton_raphson
        delta_split_phase_case(
            DELTA_PATH + "split-phase-flat-start",
            all_const_y=False,
            params=dict(
                pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-8), calculation_initialization="flat"
            ),
        ),
        delta_split_phase_case(
            DELTA_PATH + "split-phase-linear",
            all_const_y=True,
            params=pf_params(["linear", "linear_current"], 1e-8, 1e-8),
        ),
        small_delta_case(
            DELTA_PATH + "symmetric-not-supported",
            ("sym_output", "ReducedPhaseNotSupported"),
            raises_params(
                ["newton_raphson", "iterative_current"],
                "ReducedPhaseNotSupported",
                "Symmetric power flow is not supported in a network with delta-connected loads or generators.",
            ),
        ),
        Case(
            path=DELTA_PATH + "delta-on-absent-phase",
            u_rated=10e3,
            nodes=[1, 2, 3],
            lines=[Line(101, 1, 2, "abc", 2.0), Line(102, 2, 3, "ab", 3.0)],
            loads=[Load(301, 3, "asym_load", "const_pq", [0.0, 50e3, 0.0], [0.0, 10e3, 0.0], connection="delta")],
            params=raises_params(
                ["newton_raphson", "iterative_current"],
                "LoadGenOnAbsentPhase",
                "Load or generator 301 has a non-zero specified power on phases (c) that are not present at its"
                " node (node phases: ab). A delta (phase-to-phase) entry needs both of its phases: ab needs a and"
                " b, bc needs b and c, ca needs c and a.",
            ),
            raises_only=("asym_output", "LoadGenOnAbsentPhase"),
        ),
        small_delta_case(
            "state_estimation/asym-load-delta-not-supported",
            ("asym_output", "ReducedPhaseNotSupported"),
            raises_params(
                ["iterative_linear", "newton_raphson"],
                "ReducedPhaseNotSupported",
                "State estimation is not supported in a network with delta-connected loads or generators.",
            ),
        ),
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
        params=pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-8),
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
        params=pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-8),
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
        params=pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-7),
    ),
    Case(
        path="power_flow/asym-line/reduced-phase-feeder",
        u_rated=10e3,
        nodes=FEEDER_NODES,
        lines=feeder_lines(),
        loads=feeder_loads(all_const_y=False),
        params=pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-6),
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
            "calculation_method": ["newton_raphson", "iterative_current"],
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
    *swer_short_circuit_cases(),
    # the common mode of n3 is held to earth only by its loads (and the 1e-8 low susceptance); iterative_current
    # converges because, in a network with an absent phase, it puts the const_y loads into its matrix (M11c), which
    # anchors the common mode; only the small const_pq part stays explicit (11 iterations at IC_ERROR_TOLERANCE)
    phase_to_phase_case(
        "power_flow/single-phase-transformer/phase-to-phase",
        all_const_y=False,
        methods=["newton_raphson", "iterative_current"],
    ),
    phase_to_phase_case(
        "power_flow/single-phase-transformer/phase-to-phase-linear",
        all_const_y=True,
        methods=["linear", "linear_current"],
    ),
    # the network has const_pq and const_i loads, so that newton_raphson (linear start) and iterative_current (the
    # average-source start with the split-phase slot pattern) iterate; the flat start is a case of its own
    split_phase_case(
        "power_flow/single-phase-transformer/split-phase",
        all_const_y=False,
        params=pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-8),
    ),
    # iterative_current ignores calculation_initialization and always uses the average-source start, which equals the
    # flat start here (one source, u_ref 1.0); the flat start itself is exercised by newton_raphson
    split_phase_case(
        "power_flow/single-phase-transformer/split-phase-flat-start",
        all_const_y=False,
        params=dict(pf_params(["newton_raphson", "iterative_current"], 1e-8, 1e-8), calculation_initialization="flat"),
    ),
    split_phase_case(
        "power_flow/single-phase-transformer/split-phase-linear",
        all_const_y=True,
        params=pf_params(["linear", "linear_current"], 1e-8, 1e-8),
    ),
    split_phase_short_circuit_case(),
    Case(
        path="power_flow/single-phase-transformer/missing-phase-for-winding",
        u_rated=SWER_U,
        nodes=[1, 2, 3],
        lines=[Line(101, 1, 2, "a", 5.0)],
        transformers=[Transformer1P(151, 2, 3, "ab", "a", 22e3, 12.7e3, 200e3, 0.05, 2000.0)],
        params={
            "calculation_method": ["newton_raphson", "iterative_current"],
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
    *delta_cases(),
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
