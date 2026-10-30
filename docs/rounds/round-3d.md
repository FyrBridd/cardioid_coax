# Round 3d: the passive and semi-active crossovers of bass260-580

Date: 2026-10-30. The owner wanted to know what the design costs if it is run from one amplifier with a passive
crossover, or from two amplifier channels with a DSP only on the coax. This round designs both for bass260-580
(`params/round-3.toml`) and sets them against the fully active crossover of round 3c. Everything here is simulated.
The drivers are not bought, and nothing has been measured.

**The fully active design (round 3c) is unchanged:** bass to coax 150 Hz, coax LF to HF 1200 Hz, the same filters.
The passive and semi-active options are extra ways to drive the same box.

## What was run

- **Code:** branch `FyrBridd/passive-xo`, commit `3accce0`. The BEM fields came from the cache, so there were no new solves.
  Outputs are not committed.
- **Command:** `uv run speaker passive --params params/round-3.toml --out out/passive/round-3-design --trusted-limit 2851`.
  It took 83 s and exited 0. The log is `out/passive/run.log` (it still says "tens of minutes"; that text is out of date).
- **Outputs** (`out/passive/round-3-design/`): `report.html`, `network.txt`, `parts.csv`, `semi_active.txt`,
  `scores.json`, `schematic_bass.svg`, `schematic_coax.svg`, `links.svg`, `vituixcad_netlist.txt` and `plots/`.

## What changed from the plan while building it

1. **The level is tied to the woofers'.** The fit's level must stay within 1 dB of the woofers' own level. A first
   preview with the level left free threw away about 7 dB of sensitivity through a 51 Ω series resistor.
2. **The woofers must stay quiet under the coax LF.** From twice the bass crossover (300 Hz) up to 600 Hz they must
   sit at least 10 dB under the coax LF. Without this, the passive woofers ran about 8 dB under the coax out to
   about 800 Hz. That filled the coax LF's low-end droop (the cardioid chamber leaves the coax LF about 9 dB down
   at 150 Hz) and let non-cardioid sound into the band the cardioid section must own.
3. **A coax LF contour (shelf) network is offered, but not used.** It was tried first. On this design the search left
   it out: it improved flatness by only 0.08 dB.
4. **Semi-active: the coax DSP channel gets EQ.** The design said gain only. The EQ is fitted before the DSP
   high-pass, as the active fit does, and only inside 150 Hz to 2 times the LF-HF point (150 Hz to 2.4 kHz).
5. **Terminal links K1 to K5.** Each network is isolated at its boundary, so the same boards serve all three modes
   (`links.svg`, and the table at the end of `network.txt`).
6. **Impedance and tolerance rules.** The floor is 3.5 Ω. The tolerance check draws every inductor and capacitor at
   plus or minus 5 % (200 draws). The series-woofer fallback is tried only if no parallel network clears 3.5 Ω.

## The passive network

**Crossover points: bass to coax 150 Hz, coax LF to HF 1200 Hz.** The search took 30 grid pairs and kept these.
It uses the woofers in parallel, the bass inverted and the coax HF inverted. One amplifier, one pair of terminals
(IN), and links K1 to K5 all closed. Every value is a standard value (capacitors E12, resistors E24).

| Part | Value | Type | DCR | Resistor rating |
|---|---|---|---|---|
| L_b1 | 8 mH | cored inductor | 0.40 Ω | |
| C_b1 | 3.3 µF | film capacitor | | |
| C_h1 | 120 µF | film capacitor | | |
| L_h1 | 20 mH | air-core inductor | 2.24 Ω | |
| L_l1 | 3.5 mH | air-core inductor | 0.94 Ω | |
| C_l1 | 5.6 µF | film capacitor | | |
| L_ln | 0.1 mH | air-core inductor | 0.16 Ω | |
| C_ln | 33 µF | film capacitor | | |
| R_ln | 47 Ω | resistor | | 5 W |
| R_ls | 22 Ω | resistor | | 5 W |
| R_lp | 16 Ω | resistor | | 5 W |
| C_t1 | 3.3 µF | film capacitor | | |
| L_t1 | 3.5 mH | air-core inductor | 0.94 Ω | |
| C_t2 | 15 µF | film capacitor | | |
| R_ts | 0.1 Ω | resistor | | 5 W |
| R_tp | 2.4 Ω | resistor | | 10 W |

The bass is a 2nd-order low-pass (L_b1, C_b1). The coax LF has a high-pass (C_h1 with L_h1), a low-pass
(L_l1, C_l1), a notch at 2700 Hz (L_ln, C_ln, R_ln) and a pad (R_ls, R_lp). The coax HF has a high-pass
(C_t1, L_t1, C_t2) with a pad (R_ts, R_tp). The 120 µF C_h1 is a film capacitor in the coax network. Parts per pair
are in `parts.csv`.

**What the search kept and left out:**
- Kept the coax LF notch at 2700 Hz (flat23 1.55 to 1.43 dB).
- Left out the coax LF contour (1.55 to 1.47 dB, below the 0.1 dB needed), the Zobel across the woofers, a 3rd-order
  bass low-pass, a 3rd-order coax LF low-pass, and the RLC shunt across the coax HF (it made LF-HF DI smoothness
  worse by 0.37 dB, 0.25 dB allowed).

**Impedance and level:**
- Minimum impedance is **3.61 Ω at 20 kHz**, which is 0.11 Ω over the 3.5 Ω floor (`plots/impedance.png` shows the
  3.5 Ω line). The floor is the amplifier's minimum. The woofers on their own (semi-active bass channel) fall to
  3.26 Ω at 149 Hz, so that channel is below 3.5 Ω.
- Sensitivity is **79.5 dB at 2 m for 2.83 V** (23°, 300-1000 Hz). That is 1.3 dB below what the woofers give alone.
- **Tolerance check passed.** The worst of 200 draws was 3.49 Ω (limit 3.3 Ω). Flatness (flat23) was 1.43 dB nominal and
  1.65 dB at the 95th percentile (limit +0.5 dB). Whole_rear ran 10.19 to 10.67 dB (5th to 95th percentile).
- No series-woofer fallback was needed: a parallel network cleared the floor.

**Power.** Resistor ratings follow the rule in `network.txt`: the first of 5, 10, 15, 20, 30, 50 W at or above the
dissipation with one third of the amplifier's voltage (roughly what music averages). Two are close to their ratings:
R_ls at 4.90 W against 5 W, and R_tp at 9.69 W against 10 W. At full voltage they would dissipate 44 W and 87 W.
I would buy higher-rated parts than the listed ones, for margin. These are sine figures, so they are conservative for music.

### Warnings in network.txt

Two parts sit at the edge of what the search was allowed to try:
- **L_h1 is at the 20 mH limit.** The fit wanted an even larger coil for the coax LF high-pass, but 20 mH air-core is
  already 2.24 Ω of DCR and would be large and costly. The response near this part is a compromise.
- **R_ts is at the 0.1 Ω lower limit.** It barely acts. The fit wanted it out of the way, so this costs little, and
  the part could probably be a link.

## The semi-active DSP

The bass channel drives the bare woofers in parallel. The coax channel feeds the passive coax boards, with the high-pass
(C_h1, L_h1) out of circuit. A 150 Hz LR4 splits them. Use 96 kHz; 48 kHz biquads are also in `semi_active.txt`.

| Way | Gain | Delay | Polarity | PEQ (f0, gain, Q) |
|---|---|---|---|---|
| bass | −1.78 dB | 1.35 ms | inverted | 80 Hz +1.18 dB Q1.65; 214 Hz −1.77 dB Q1.22 |
| coax | 0 dB | 0 | normal | 153 Hz +8.54 dB Q1.36; 179 Hz −0.84 dB Q6.0 (limit); 752 Hz −0.80 dB Q0.5 (limit) |

- **Q limits.** Two coax PEQs sit at the fitter's Q limits: the 179 Hz one at the top (6.0, a very narrow cut) and the
  752 Hz one at the bottom (0.5, very broad). Both are small (0.8 dB). The marker in `semi_active.txt` only checks
  gain, so it does not flag them. The fit is working at the edge of what the PEQs can do, so these two do little.
- **Stacked boost.** The coax channel's combined boost peaks at 8.2 dB at 151 Hz. That is large enough to need the
  usual care. Amplifier clipping and the coax LF's excursion at 150 Hz should be checked (see Open). The bass gain is
  −1.78 dB against the coax. Normalise in the DSP so the largest boost does not clip, as the file says.

## Three modes, trusted scores

Trusted band (150 Hz to 2851 Hz; bass DI from 120 Hz). Lower is better except the rear and cardioid rows.

| Score (dB) | Passive | Semi-active | Active |
|---|---|---|---|
| 23° flatness | 1.43 | 1.90 | **1.16** |
| Window spread | 3.60 | 3.50 | **3.40** |
| Presence lift, 2-2.85 kHz | 4.22 | 4.33 | **3.90** |
| DI smoothness, bass, 120-500 Hz | **0.94** | 1.55 | 1.66 |
| DI smoothness, LF-HF, 500-2851 Hz | **0.77** | 0.90 | 1.08 |
| Whole-system rear rejection, 150-300 Hz (higher is better) | **10.44** | 9.96 | 9.88 |
| Cardioid hard / aim (coax section) | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 |

The coax section's cardioid scores are the same because the coax is the same box.

## Where passive falls short, and why

- **Flatness.** Passive is 0.27 dB behind active at 23°, and semi-active is 0.74 dB behind. I did not find out why
  semi-active is behind passive here. Both are within 0.8 dB of active. The owner can judge it from `plots/crossover_23_*.png`.
- **The 2 to 5 kHz window.** Presence lift is 4.2 dB against 3.9 dB for active. All three are well above the
  0.5 dB the score calls a pass. Passive is not the cause: the lift is the coax's, and it was already there in
  round 3c. Window spread is 3.6 dB against 3.4 dB. Above 2851 Hz the scores are low confidence and are not ranked.
- **Rear rejection looks better than it is.** The whole_rear score is the minimum over 150 to 300 Hz, so it
  is set by the 150 Hz edge, where all three are within 0.6 dB of each other. Across the band the picture is different:
  at about 250 Hz the passive rejects about 13.5 dB and semi-active and active about 18.5 dB
  (`plots/whole_rear_compare.png`). The passive's slow filters give less rear rejection in the middle of the band.
  The score rewards passive for the one frequency where all three are weakest.
- **Sensitivity.** 79.5 dB is 1.3 dB below the woofers alone. It needs about 1.3 dB more amplifier drive than the woofers alone for the same level. The active figure is not comparable, as each way has its own amplifier gain.
- **Size and cost.** A 20 mH air-core coil, a 120 µF film capacitor and two 87 W-at-full-voltage resistors are not small
  or cheap. Parts cost was not estimated.

## Notes

- **Models.** The parts values come from the model of unmeasured drivers. Real driver impedance will move the
  network. Treat it as a first draft to be trimmed after measurement.
- **Above 2851 Hz** the coax HF is low confidence. The coax HF network was fitted on the trusted band only.
- **Polarity.** In passive mode the bass and coax HF are wired inverted (+ terminal to the network's ground side). In semi-active the bass is inverted in the DSP.

## Open

- Whether the owner wants passive at all: it needs one amplifier and no DSP, at lower rear rejection through the cardioid
  band and the parts above. Semi-active needs two channels and a DSP on the coax only.
- Buy the resistors with margin, or check them against the full-voltage figures in `network.txt`.
- The stacked 8.2 dB coax boost at 151 Hz needs a check of the coax LF's excursion and amplifier headroom.
- Everything above depends on unmeasured drivers. Rerun after the coax and woofers are measured.
- The VituixCAD netlist (`vituixcad_netlist.txt`) has not been imported or checked by the owner.
- The owner's report page gets its Passive option from the controller after this round.
