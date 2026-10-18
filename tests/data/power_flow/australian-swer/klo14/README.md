<!--
SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
SPDX-FileCopyrightText: 2024 Team Nando - Public Repositories

SPDX-License-Identifier: MPL-2.0 AND BSD-3-Clause
-->

# Validation case: power flow on the Australian rural SWER network KLO14

Asymmetric power flow (`newton_raphson`) on the medium-voltage part of the rural 22 kV feeder KLO14 of AusNet Services
(Victoria, Australia), with its SWER (single-wire earth return) networks, with OpenDSS as the reference solver. It is a
batch of two scenarios: scenario 0 at minimum and scenario 1 at maximum demand. `input.json` holds the minimum demand
and `update_batch.json` sets the `asym_load` powers of both scenarios. There is no single-calculation output
(`asym_output.json`).

## Source data

- Repository: https://github.com/Team-Nando/MV-LV-Networks, commit `8429f8a19f7b95acca2965cc48f2736880e5b31c`
- File: `Network_2_Rural_KLO14.xlsx`, sha256 `19d745c6b6aeedf78b2142d5158b4d7761e9e8b0dd35f4f311156a2187cdc46a`
- Licence: BSD 3-Clause, Copyright (c) 2024, Team Nando - Public Repositories (`LICENSES/BSD-3-Clause.txt`). The files of this case are
  derived from it and are licensed `MPL-2.0 AND BSD-3-Clause` (see the `.license` files).
- The model is derived from the repository's `backend.py` (its OpenDSS model of the network), with the deviations
  listed below. The workbook is not part of this repository.

## Network

| Element | Count |
|---|---|
| nodes | 3224: the 66 kV source bus, 2156 22 kV, 367 SWER and 700 LV busbars |
| `line` (three-phase, 22 kV) | 2154 |
| `asym_line` (one-phase SWER, phase a) | 360 |
| zone substation (`transformer`) | 66/22.7 kV, 100 MVA, Dyn1 |
| voltage regulators (`transformer`, YNyn0) | 3: PLENTY_VLY_RISE (200 kVA), MERRIANG_JANNA_REG (5000 kVA), DONNYBROOK_RD_REG (10000 kVA) |
| ISO transformers (`single_phase_transformer`) | 7: EDEN_PARK_STH_ISO (100 kVA), TWELTH_AVE_ISO (100 kVA), EDEN_PARK_NORTH_ISO (100 kVA), NINTH_AVE_ISO (100 kVA), UPPER_PLENTY_STH_ISO (200 kVA), UPPER_PLENTY_ISO (100 kVA), FOURTEENTH_AVE_ISO (100 kVA) |
| three-phase distribution transformers (`transformer`, 22/0.433 kV Dyn1) | 137 |
| split-phase distribution transformers (`single_phase_transformer`, 22 kV / 250-0-250 V) | 358 |
| SWER distribution transformers (`single_phase_transformer`, 12.7/0.24 kV) | 161 |
| 22 kV phase-to-earth distribution transformers (`single_phase_transformer`, 12.7/0.25 kV) | 44 |
| capacitors (`shunt`, 900 kvar) | 2 |
| loads (`asym_load`, const_pq) | 700 |

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
  phase-to-phase single-phase winding. OpenDSS gets matching shunt Reactors (`Reactor.lo{id}n{node}`), whose powers
  are added to the from side of the transformer in the reference.
- The source has Z1 = 0.52824 + j2.113 ohm and Z0 with the magnitude of backend's 0.59157 + j1.7 ohm but Z1's angle
  (PGM has a real `z01_ratio`). This is irrelevant: the zone substation primary is delta, so the 66 kV source carries
  no zero-sequence current.
- Line rows named `Delete` are dropped. One-phase SWER lines use Z1 and C1 of their LineCode only, as OpenDSS solves
  backend's 1-phase LineCodes (r0, x0 and b0 are ignored), so they have no earth-return reactance in either tool. All
  22 kV and 12.7 kV SWER nodes have `u_rated` 22 kV (the SWER convention of power-grid-model).
- Capacitors are grounded-wye shunts with b1 = b0 = kvar / kV^2.
- Near-zero lines are merged (in both tools): see below.
- OpenDSS runs with tolerance 1e-11 (1e-12 stalls at round-off on SMR8) and MaxIterations
  200.

## Merged near-zero lines

The source data has connector lines with a series impedance down to 1.4e-7 ohm and a zone substation with uk =
1.41e-6. They make the admittance matrix so ill-conditioned that PGM's Newton-Raphson stalls at a max deviation of
about 3e-7 (error tolerance 1e-8) and PGM and OpenDSS agree to only about 2e-7. This is not caused by the reduced-phase
or SWER parts: the three-phase part of the network on its own (the upstream code path) stalls at the same level.

Therefore every line whose series impedance |z1| L is below 0.001 ohm (2.1e-06 pu on the 22 kV base of 484 ohm,
1 MVA) is removed and its two end nodes are merged into one, in both tools. A merged group of nodes keeps the node
closest to the source (the fewest branches; on a tie the lowest id), and every element on a merged node moves to it:
lines, transformer terminals and capacitors. A three-phase line merges 22 kV nodes and a one-phase line SWER nodes;
the generator asserts that a 22 kV node never merges with a SWER node. A kept line whose ends merge is dropped (it runs
in parallel with merged lines only), and lines that end up between the same two nodes stay parallel. Every other line
keeps the length and LineCode of the data. The zone substation uk is raised to 0.0001 (uk and pk scaled together,
r / x kept): with a smaller uk it is the smallest impedance left after the merge, and PGM misses the tolerances on
KLO14. Here: 1105 three-phase and 205 one-phase lines are merged (|z1| L up to 0.000997 ohm), and 1309 nodes with them, which leaves 3224 of 4533 nodes and 2154 of 3259 three-phase and 360 of 565 one-phase lines; 0 kept lines became a self-loop and were dropped, and 0 lines run in parallel after the merge. The smallest |z1| L left is 0.00109 ohm. The zone substation uk is raised from 1.41e-06 to 0.0001 (0.00052 ohm on the 22.7 kV side).

The merge changes the OpenDSS solution against the source network (every line of the data and the sheet uk of the zone
substation, solved with the same tolerance) by the following. A node voltage change is the largest
|V - V_source| / |V_source| over the phases of the kept nodes of a class (the phasors, so the angle counts too); a
merged node is compared at the node it is merged into.

| change against the source network | scenario 0: minimum demand | scenario 1: maximum demand |
|---|---|---|
| node voltage, 66 kV source | 1.3e-06 | 1.2e-06 |
| node voltage, 22 kV | 5.5e-05 | 5.2e-05 |
| node voltage, SWER 12.7 kV | 4.9e-05 | 4.3e-05 |
| node voltage, LV 433 V three-phase | 5.5e-05 | 5.1e-05 |
| node voltage, LV 250-0-250 V split-phase | 5.2e-05 | 4.9e-05 |
| node voltage, LV 240 V SWER | 4.9e-05 | 4.3e-05 |
| node voltage, LV 250 V 22 kV phase-to-earth | 3.9e-05 | 3.2e-05 |
| node voltage, merged nodes | 5.5e-05 | 5.2e-05 |
| losses | 91.851 to 91.762 kW (-9.7e-04) | 303.982 to 303.865 kW (-3.8e-04) |
| feeder P (relative) | -2.3e-05 | -1.6e-05 |

## Demand

One const_pq `asym_load` per distribution transformer at power factor 0.95 lagging, a fraction of its kVA:
0.1 at minimum and 0.2 at maximum demand. Three-phase loads are split equally over the phases,
split-phase loads over the two legs (to the centre tap), and phase-to-earth loads are on phase a. At maximum demand the
SWER loads behind an ISO get the coincidence factor min(1, S_iso / (0.2 x the kVA of the SWER transformers
behind it)); applied here: none (no ISO is above its rating). A higher maximum (the 0.60 of the original plan) does not converge in OpenDSS:
the voltage collapses on both networks.

## OpenDSS results

| | scenario 0: minimum demand (0.1) | scenario 1: maximum demand (0.2) |
|---|---|---|
| OpenDSS iterations | 9 | 15 |
| feeder P, Q | 3.460 MW, -1.644 Mvar | 7.041 MW, 0.031 Mvar |
| load P | 3.369 MW | 6.737 MW |
| losses | 91.8 kW | 303.9 kW |
| 66 kV source (pu of 38105.1 V) | 1.0002 to 1.0006 | 0.9988 to 0.9995 |
| 22 kV (pu of 12701.7 V) | 1.0097 to 1.0429 | 0.9483 to 1.0312 |
| SWER 12.7 kV (pu of 12700.0 V) | 1.0041 to 1.0383 | 0.9333 to 1.0114 |
| LV 433 V three-phase (pu of 250.0 V) | 0.9804 to 1.0362 | 0.9285 to 1.0207 |
| LV 250-0-250 V split-phase (pu of 250.0 V) | 0.9790 to 1.0313 | 0.9101 to 1.0225 |
| LV 240 V SWER (pu of 240.0 V) | 1.0241 to 1.0622 | 0.9512 to 1.0316 |
| LV 250 V 22 kV phase-to-earth (pu of 250.0 V) | 1.0176 to 1.0192 | 0.9574 to 0.9608 |
| highest ISO loading | 46.4% (UPPER_PLENTY_ISO) | 95.3% (UPPER_PLENTY_ISO) |
| highest regulator loading | 46.8% (PLENTY_VLY_RISE) | 78.0% (PLENTY_VLY_RISE) |

## Reference output

`asym_output_batch.json` has, per scenario, `node` (u, u_angle; 0 on absent phases), `asym_line` (p and q, from and
to), `single_phase_transformer` and `transformer` (p and q, from and to, and the loading), rounded to 12 significant
digits. It has no three-phase `line` output: the node voltages at both ends pin every line, and it keeps the data small.
All datasets are written in the compact JSON format (a root `attributes` list per component and list rows).

Tolerances: rtol 1e-06 and atol 1e-08 (0.001 on `[pq]_(from|to)`), the smallest round values with a 10x margin over both scenarios, with the
atol of the power fields at most 1e-3 (W, var). The power fields have their own atol: a branch power is computed from
the voltages at its ends, so a short branch amplifies their round-off. Measured: newton_raphson converges in 3 and 4 iterations (minimum, maximum demand); the node voltages agree to 9.8e-11 relative and the powers to 1.2e-3 W absolute, a worst ratio of 0.069. iterative_current (7 and 9 iterations) misses the rule (worst ratio 0.207, on the powers through the regulator PLENTY_VLY_RISE), so it is not run.

## IDs

An MV node keeps its AusNet Node_ID; the 66 kV source bus is 1000000; the LV busbar of distribution transformer E is
2000000 + E; a line, transformer, capacitor or the zone substation (Element_ID E) is 3000000 + E; the load of
distribution transformer E is 4000000 + E; the source is 5000000. In OpenDSS every bus is `n{node id}`, and the
elements are `Line.l{id}`, `Transformer.t{id}`, `Capacitor.c{id}` and `Load.ld{id}_{a,b,c}`.

## Generation

The files are generated by `tests/data/generators/australian_swer_networks.py` with OpenDSS (opendssdirect.py
0.9.4) as the reference solver; do not edit them by hand. `test.dss` is exactly the script that was run:
the generator runs it line by line and reads each reference right after the `Solve` of its scenario. From the
repository root (the workbooks are downloaded from the pinned commit, or read from a local directory with
`--source-dir`):

```sh
uv run --no-project --with opendssdirect.py --with numpy --with openpyxl \
    python tests/data/generators/australian_swer_networks.py tests/data klo14
```
