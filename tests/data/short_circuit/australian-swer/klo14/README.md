<!--
SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
SPDX-FileCopyrightText: 2024 Team Nando - Public Repositories

SPDX-License-Identifier: MPL-2.0 AND BSD-3-Clause
-->

# Validation case: short circuit on the Australian rural SWER network KLO14

IEC 60909 short circuit (maximum, c = 1.1) on the medium-voltage part of the rural 22 kV feeder KLO14 of AusNet
Services (Victoria, Australia), with its SWER (single-wire earth return) networks, with OpenDSS as the reference
solver. The network is the one of the power-flow case `tests/data/power_flow/australian-swer/klo14`, without loads.
It is a batch of 9 fault scenarios: `input.json` has the fault (one ID, 6000000) of scenario 0, and
`update_batch.json` points it to the node of each scenario (fault_object, fault_type, fault_phase, r_f, x_f). There is
no single-calculation output (`sc_output.json`).

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

## Modelling and deviations from backend.py

- The LV networks are replaced by one lumped load per distribution transformer on a new LV busbar. The medium-voltage
  network, the SWER networks and every transformer of the source data are kept.
- Voltage regulators are fixed at the neutral tap as a YNyn0 bank on the rating with backend's xhl and %loadloss.
  backend builds each one as three 12.7/1.27 kV autotransformers under a RegControl starting at tap 0; at the neutral
  tap they are a 1:1 connection plus the leakage impedance, which both tools model identically as YNyn0 (in positive
  and zero sequence).
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

The faults are selected on the merged network, so every fault node exists in both tools (a merged node does not). The
merge changes the OpenDSS fault currents against the source network (every line of the data and the sheet uk of the
zone substation, with the same faults) by the following: the largest |I_f - I_f,source| over the phases, relative to
the largest |I_f,source|.

| scenario | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| change of I_f | 1.8e-03 | 2.3e-05 | 3.5e-04 | 2.7e-05 | 2.4e-05 | 3.2e-05 | 2.7e-05 | 1.6e-05 | 2.0e-05 |

## Fault scenarios

One fault (ID 6000000) per scenario, selected deterministically by the generator. The path impedance of a node is
the sum |z1| L of the lines on its path from the source bus along the breadth-first tree of the network (the branches
of a node taken in the order of their IDs; transformers count 0); every tie goes to the lowest ID. A bolted fault
(r_f = 0) is a Fault of 1e-13 ohm in OpenDSS (OpenDSS replaces r = 0 by 1e-4 ohm), a fault with x_f a series R + jX
Reactor.

| # | node | where | fault | r_f, x_f (ohm) | why this node | OpenDSS \|I_f\| (A) |
|---|---|---|---|---|---|---|
| 0 | 3849 | zone-substation 22 kV bus | three_phase abc | 0, 0 | the secondary of the zone substation | 50788, 50788, 50788 |
| 1 | 2088 | 22 kV line end | single_phase_to_ground a | 0, 0 | the 22 kV node with the largest path impedance from the zone-substation bus (91.92 ohm) | 163.516 |
| 2 | 741 | secondary bus of the regulator DONNYBROOK_RD_REG | two_phase bc | 0.5, 0 | the largest regulator (10000 kVA); r_f = 0.5 ohm, since a bolted 22 kV bridge has a round-off of about 1.6e-2 in OpenDSS | 9176.24, 9176.24 |
| 3 | 3044 | SWER line end | single_phase_to_ground a | 0, 0 | the SWER node with the largest path impedance from its ISO secondary (UPPER_PLENTY_STH_ISO, 18.89 ohm) | 269.216 |
| 4 | 3044 | SWER line end (high-impedance fault, a Reactor in OpenDSS) | single_phase_to_ground a | 20, 5 | the node of scenario 3 | 196.863 |
| 5 | 2003855 | LV busbar of a SWER transformer (240 V) | single_phase_to_ground a | 0, 0 | the SWER transformer 3003855 (10 kVA) whose primary has the largest path impedance from its ISO secondary (UPPER_PLENTY_STH_ISO, 18.89 ohm) | 1379.05 |
| 6 | 2004222 | LV busbar of a split-phase transformer (250-0-250 V) | two_phase ab | 0, 0 | the largest split-phase transformer, 3004222 (100 kVA); a bolted two-phase fault is accurate in OpenDSS here (the shorted legs sit near 0 V) | 5035.53, 5035.53 |
| 7 | 2004290 | LV busbar of a three-phase transformer (433 V) | three_phase abc | 0, 0 | the largest three-phase transformer, 3004290 (750 kVA) | 16080.6, 16071, 16076.3 |
| 8 | 2004010 | LV busbar of a 22 kV phase-to-earth transformer (250 V) | single_phase_to_ground a | 0, 0 | the largest 22 kV phase-to-earth transformer, 3004010 (25 kVA) | 2939.92 |

Self-checks of the generator: OpenDSS converges in every scenario; the Reactor of the high-impedance fault carries
V / z_f with the OpenDSS node voltage (scenario 4: 1.5e-16 relative); the bolted two-phase fault, solved again on its own with a 10x
larger bridge, gives the same fault current (scenario 6: 9.9e-12 relative).

## Short-circuit model

- The source is c times 66 kV behind c times its impedance (IEC 60909, as in PGM), with c = 1.1 for the maximum
  (PGM takes it from the u_rated of the source node).
- No loads. The capacitors and the line capacitances stay in, as in PGM's short-circuit admittance matrix.
- `test.dss` builds the network once; per scenario it defines the fault elements (`Fault.f6000000s{k}` or
  `Reactor.f6000000s{k}` for scenario k), solves, and disables them again. OpenDSS runs with tolerance
  1e-11 in snapshot mode.

## Reference output

`sc_output_batch.json` has, per scenario, `fault` (i_f, i_f_angle), `single_phase_transformer` and `transformer` (i
and i_angle, from and to; the current of a low-susceptance Reactor is added to the from side), rounded to 12
significant digits. It has no `node`, `line` or `asym_line` output: the transformer currents pin the network around each fault, and the data of the four cases on both networks has a size budget of 10 MB (it is 9.0 MB; the asym_line output would add 1.05 MB and the node output 5 MB). All datasets are written in the
compact JSON format (a root `attributes` list per component and list rows).

Tolerances: rtol 1e-08 and atol 1e-08 (1e-06 on `i_(from|to)(_angle)?`), the smallest round values with a 10x margin over all scenarios. Measured: the fault currents agree to 9.2e-11 relative (as phasors, so the angles too), the branch currents above 1 A to 6.0e-9 relative and those up to 10 A to 3.2e-8 A absolute, a worst ratio of 0.030 (the 66 kV current of the zone substation in scenario 3, 3.8 A).

## IDs

As in the power-flow case: an MV node keeps its AusNet Node_ID; the 66 kV source bus is 1000000; the LV busbar of
distribution transformer E is 2000000 + E; a line, transformer, capacitor or the zone substation (Element_ID E) is
3000000 + E; the source is 5000000 and the fault 6000000. In OpenDSS every bus is `n{node id}`, and the elements
are `Line.l{id}`, `Transformer.t{id}` and `Capacitor.c{id}`.

## Generation

The files are generated by `tests/data/generators/australian_swer_networks.py` with OpenDSS (opendssdirect.py
0.9.4) as the reference solver, together with the power-flow case; do not edit them by hand. `test.dss` is
exactly the script that was run: the generator runs it line by line and reads each reference right after the `Solve`
of its scenario. From the repository root (the workbooks are downloaded from the pinned commit, or read from a local
directory with `--source-dir`):

```sh
uv run --no-project --with opendssdirect.py --with numpy --with openpyxl \
    python tests/data/generators/australian_swer_networks.py tests/data klo14
```
