# Round 3d: the passive and semi-active crossovers of bass260-580

Date: 2026-10-30. The owner wanted to know what the design costs if it is run from one amplifier with a passive
crossover, or from two amplifier channels with a DSP only on the coax. This round designs both for bass260-580
(`params/round-3.toml`) and sets them against the fully active crossover of round 3c. Everything here is simulated.
The drivers are not bought, and nothing has been measured.

**The fully active design (round 3c) is unchanged:** bass to coax 150 Hz, coax LF to HF 1200 Hz, the same filters.
The passive and semi-active options are extra ways to drive the same box.

**Revised after the whole-branch review.** The first version of this write-up had a bass network that resonated with
the woofers' voice coils at about 2.4 kHz, where the model could not see it (item 7 below). The fit now has a rule
against that, the links were made two-pole, and every number below is from the rerun.

## What was run

- **Code:** branch `FyrBridd/passive-xo`, commit `df1bd75`. The BEM fields came from the cache, so there were no new
  solves. Outputs are not committed.
- **Command:** `uv run speaker passive --params params/round-3.toml --out out/passive/round-3-design --trusted-limit 2851`.
  It took 129 s and exited 0. The log is `out/passive/run.log`. `--trusted-limit` is now required: pass the trusted
  limit of the directivity run.
- **Outputs** (`out/passive/round-3-design/`): `report.html`, `network.txt`, `parts.csv`, `semi_active.txt`,
  `scores.json`, `schematic_bass.svg`, `schematic_coax.svg`, `links.svg`, `vituixcad_netlist.txt` and `plots/`
  (now with `electrical.png`).

## What changed from the plan while building it

1. **The level is tied to the woofers'.** This is a soft target, not a hard rule. The fit's level is the weighted
   mean of the 23° response (1/3-octave smoothed) over 100 Hz to 10 kHz: weight 1 from 100 Hz to the trusted limit
   (2851 Hz), and 0.25 from there to 10 kHz. The fit is penalised when that level falls more than 1 dB under the
   woofers' own level (80-150 Hz). **In this run it sits 3.5 dB under**, so the anchor gave way by 2.5 dB: the
   stopband rule (item 7) cost level. The sensitivity is measured differently (23°, 300-1000 Hz), and is 3.3 dB under
   what the woofers give alone. A first preview run (not in this round's outputs) with the level left free threw away
   about 7 dB through a 51 Ω series resistor.
2. **The woofers must stay quiet under the coax LF.** From twice the bass crossover (340 Hz) up to 600 Hz they
   should sit at least 10 dB under the coax LF. This is a penalty in the fit, not a guarantee. Without it, the passive
   woofers ran about 8 dB under the coax out to about 800 Hz (an earlier preview run). That filled the coax LF's
   low-end droop (the cardioid chamber leaves the coax LF about 9 dB down at 150 Hz) and let non-cardioid sound into
   the band the cardioid section must own. **The delivered network nearly meets it, and misses at the bottom of the
   band.** Recomputed from this run (23°, 1/3-octave smoothed, as the fit does): the woofers sit 8.8 dB under the
   coax LF at 356 Hz and 9.7 dB at 378 Hz, so 1.2 dB short at worst, 2 of 11 grid points. From 400 Hz up they are
   10.8 dB under or more, reaching 21.4 dB at 600 Hz.
3. **A coax LF contour (shelf) network is offered, and this time it is used.** It is R_lc (5.6 Ω) with L_lc
   (0.7 mH) across it, before the coax LF pad. It improved flatness from 1.91 to 1.75 dB.
4. **Semi-active: the coax DSP channel gets EQ.** The design said gain only. The EQ is fitted before the DSP
   high-pass, as the active fit does, and only inside 150 Hz to 2 times the LF-HF point (150 Hz to 2.4 kHz).
5. **Terminal links K1 to K5.** Each network is isolated at its boundary, so the same boards serve all three modes
   (`links.svg`, and the table at the end of `network.txt`). **K3, K4 and K5 are two-pole:** each breaks both the hot
   and the return between a network and its driver. In the active modes every terminal pair then reaches its driver
   alone. Without that, the drivers' returns stayed joined at the network ground, and in passive mode the inverted
   woofers have their + on that ground: wiring the BASS terminals "straight to the woofers" could have shorted an
   amplifier channel.
6. **Impedance and tolerance rules.** The floor is 3.5 Ω. The tolerance check draws every inductor and capacitor at
   plus or minus 5 % (200 draws). The series-woofer fallback is tried only if no parallel network clears 3.5 Ω.
7. **New: an electrical stopband rule.** The first network (L_b1 8 mH, C_b1 3.3 µF) had a bass transfer of
   **+3.0 dB at 2.43 kHz**: C_b1 was far too small to shunt the woofers, and it resonated with their voice-coil
   inductance. The fit could not see it. The bass field is exactly zero above the trusted limit, the bass-under term
   stops at 600 Hz, and the flatness is 1/3-octave smoothed. Now each branch's voltage transfer, on the modelled
   driver loads and independent of the acoustic data, must stay 12 dB under its level at its crossover point, from
   4 times past that point on: the bass from 4 × 170 = 680 Hz up, the coax LF from 4.8 kHz up and below 42 Hz, and
   the coax HF below 300 Hz. The fit penalises every dB over that at 20 per dB (`STOP_WEIGHT`). On this design weights
   3 and 10 left the bass 0.7 and 1.8 dB over the line near 4 times its crossover point. Weights 15 to 100 held it
   within 0.3 dB. 20 is one step above the lowest that held. `plots/electrical.png` shows each branch with the
   limits.

   **What the rule costs here.** Flatness at 23° goes from 1.43 dB (first network) to 1.75 dB. The bass-to-coax
   point moves from 150 Hz to 170 Hz. The level drops by about 2 dB (item 1). It does not force the series fallback.

## The passive network

**Crossover points: bass to coax 170 Hz, coax LF to HF 1200 Hz.** The search took 30 grid pairs and kept these.
It uses the woofers in parallel, the bass inverted and the coax HF inverted. One amplifier, one pair of terminals
(IN), and links K1 to K5 all closed. Every value is a standard value (capacitors E12, resistors E24).

| Part | Value | Type | DCR | Peak current at full voltage | Resistor rating |
|---|---|---|---|---|---|
| L_b1 | 15 mH | cored inductor | 0.54 Ω | 9.6 A (saturation current at least this) | |
| C_b1 | 100 µF | bipolar electrolytic capacitor | | | |
| C_h1 | 120 µF | film capacitor | | | |
| L_h1 | 20 mH | air-core inductor | 2.24 Ω | 6.7 A | |
| L_l1 | 2.5 mH | air-core inductor | 0.79 Ω | 7.3 A | |
| C_l1 | 10 µF | film capacitor | | | |
| R_lc | 5.6 Ω | resistor | | | 5 W |
| L_lc | 0.7 mH | air-core inductor | 0.42 Ω | 1.4 A | |
| R_ls | 47 Ω | resistor | | | 10 W |
| R_lp | 56 Ω | resistor | | | 5 W |
| C_t1 | 3.3 µF | film capacitor | | | |
| L_t1 | 1 mH | air-core inductor | 0.50 Ω | 2.3 A | |
| C_t2 | 22 µF | film capacitor | | | |
| R_ts | 0.82 Ω | resistor | | | 5 W |
| R_tp | 1.6 Ω | resistor | | | 10 W |

The peak current is the peak of a sine with the amplifier's full 25 V rms on the input, the largest from 20 Hz to
20 kHz in either mode. A cored coil must not saturate below it. Parts per pair are in `parts.csv` (new column
`peak_a`).

**What each part does.**
- **Bass (L_b1, C_b1).** A series coil and a shunt capacitor, but on the woofers' impedance it is not a textbook
  2nd-order low-pass at 170 Hz. Its transfer (`plots/electrical.png`) peaks at **+4.7 dB at 54 Hz** (the coil and
  capacitor resonating, with the woofers' sealed-box impedance peak near there). It is −13.7 dB at 170 Hz and stays
  near −15 dB to about 450 Hz, a shelf rather than a roll-off. Then it falls steeply:
  −25.5 dB at 680 Hz, −34 dB at 1 kHz, and no resonance above. On the fit's frequency grid it is exactly 12.0 dB
  under its 170 Hz level at 686 Hz. On a finer grid it is 11.8 dB under at exactly 680 Hz, 0.2 dB short of the rule.
- **Coax high-pass (C_h1 from IN to node cx, L_h1 from cx to ground).** It feeds both coax networks, the LF and the
  HF. It is resonant: the coax LF gets its most voltage, −5.2 dB, at 99 Hz.
- **Coax LF.** A low-pass (L_l1, C_l1), the contour (R_lc with L_lc across it), and a pad (R_ls, R_lp).
- **Coax HF.** A high-pass (C_t1, L_t1, C_t2) and a pad (R_ts, R_tp).

**What the search kept and left out:**
- Kept the coax LF contour (flat23 1.91 to 1.75 dB).
- Left out the Zobel across the woofers, a 3rd-order bass low-pass, a 3rd-order coax LF low-pass, the coax LF notch
  at 2700 Hz and the RLC shunt across the coax HF. Each made flat23 worse.

**Impedance and level:**
- Minimum impedance is **3.52 Ω at 57 Hz**, which is 0.02 Ω over the 3.5 Ω floor. That is also the minimum in the
  bass band (20-500 Hz). `plots/impedance.png` shows the 3.5 Ω line. The floor is the SA550's minimum.
- **Semi-active impedance.** The coax channel's minimum is 3.52 Ω at 20 kHz, 0.02 Ω over the floor. The bass channel
  (bare woofers in parallel) dips to **3.26 Ω at 149 Hz, below the 3.5 Ω the SA550 needs**. The fit and the tolerance
  check do not cover that channel, because the amplifier for it is not chosen. Pick one that is happy at about 3 Ω.
- Sensitivity is **77.5 dB at 2 m for 2.83 V** (23°, 300-1000 Hz). That is 3.3 dB below what the woofers give alone
  (it was 79.5 dB, 1.3 dB below, before the stopband rule).
- **Tolerance check passed.** The worst of 200 draws was 3.43 Ω (limit 3.3 Ω). Flatness (flat23) was 1.75 dB nominal
  and 2.02 dB at the 95th percentile (limit +0.5 dB). Whole_rear ran 11.34 to 11.94 dB (5th to 95th percentile).
- No series-woofer fallback was needed: a parallel network cleared the floor.

**Power.** Resistor ratings follow the rule in `network.txt`: the first of 5, 10, 15, 20, 30, 50 W at or above the
dissipation with one third of the amplifier's voltage (roughly what music averages). Three are close to their
ratings in passive mode: R_ls at 9.52 W against 10 W, R_tp at 7.23 W against 10 W, and R_ts at 4.55 W against 5 W.
At full voltage they would dissipate 86 W, 65 W and 41 W. In semi-active mode (the same boards, but the amplifier
channel carries only the coax) R_ls reaches 9.38 W, R_tp 7.44 W and R_ts 4.68 W against 5 W. I would buy
higher-rated parts than the listed ones, at least for R_ls and R_ts. These are sine figures, so they are conservative
for music.

**Coils.** L_b1 is 15 mH on a core and must carry 9.6 A peak without saturating. That is a large, heavy part. L_h1
(20 mH air core, 2.24 Ω) is also large.

### Warnings in network.txt

- **L_h1 is at the 20 mH limit.** The fit wanted an even larger coil for the coax high-pass, but 20 mH air-core is
  already 2.24 Ω of DCR and would be large and costly. The response near this part is a compromise.

## The semi-active DSP

The bass channel drives the bare woofers in parallel. The coax channel feeds the passive coax boards (COAX + to node
cx, COAX − to node 0, the coax boards' ground), with links K1 and K2 open, so the coax high-pass (C_h1, L_h1) is out
of circuit, K3 open, and K4 and K5 closed. BASS + goes to the woofers' + terminal and BASS − to their − terminal.
A 150 Hz LR4 splits them. Use 96 kHz; 48 kHz biquads are also in `semi_active.txt`.

**Polarity** is set in the DSP, relative to each driver's own + terminal. The bass at −1 means the DSP inverts the
signal to the woofers' + terminal. The coax boards keep their own wiring through K4 and K5 (coax HF inverted, as in
passive mode). `semi_active.txt` says this in its header.

| Way | Gain | Delay | Polarity | PEQ (f0, gain, Q) |
|---|---|---|---|---|
| bass | −3.84 dB | 1.65 ms | inverted | 80 Hz +1.18 dB Q1.65; 214 Hz −1.77 dB Q1.22 |
| coax | 0 dB | 0 | normal | 150 Hz +12.0 dB Q0.72 (gain limit); 202 Hz −3.83 dB Q0.91; 610 Hz −3.43 dB Q1.94; 1991 Hz −1.40 dB Q2.53 |

- **Gain limit.** The 150 Hz coax PEQ is at the fitter's +12 dB limit, and `semi_active.txt` marks it. The new
  coax boards (without the coax high-pass) pass about 3 dB less to the coax LF at 150 to 300 Hz than the first ones
  did (−14.4 against −11.4 dB at 150 Hz), so the DSP has to fill more of the cardioid chamber's low-end droop.
- **Stacked boost.** The four coax PEQs together (my sum of their curves) peak at **9.2 dB at 139 Hz**. That needs
  the usual care: amplifier clipping and the coax LF's excursion near 140 Hz should be checked (see Open). Normalise in
  the DSP so the largest boost does not clip, as the file says.

## Three modes, trusted scores

Trusted band (150 Hz to 2851 Hz; bass DI from 120 Hz). Lower is better except the rear and cardioid rows. The first
column is the passive network before the stopband rule, for comparison.

| Score (dB) | Passive, first network | Passive | Semi-active | Active |
|---|---|---|---|---|
| 23° flatness | 1.43 | 1.75 | 2.05 | **1.16** |
| Window spread | 3.60 | 3.73 | 3.73 | **3.40** |
| Presence lift, 2-2.85 kHz | 4.22 | 4.04 | 4.63 | **3.90** |
| DI smoothness, bass, 120-500 Hz | 0.94 | **0.69** | 1.39 | 1.66 |
| DI smoothness, LF-HF, 500-2851 Hz | 0.77 | 1.11 | 1.23 | **1.08** |
| Whole-system rear rejection, 150-300 Hz (higher is better) | 10.44 | **11.62** | 10.01 | 9.88 |
| Cardioid hard / aim (coax section) | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 |

The coax section's cardioid scores are the same because the coax is the same box.

## Where passive falls short, and why

- **Flatness.** Passive is 0.59 dB behind active at 23°, and semi-active is 0.88 dB behind. The stopband rule cost
  passive 0.32 dB of it. Semi-active uses the new coax boards, whose heavier coax LF pad leaves the DSP more droop to
  fill (its 150 Hz PEQ is at its limit). The owner can judge both from `plots/crossover_23_*.png`.
- **Below the fitted band.** The fit starts at 100 Hz. Below it the passive sum dips to about 72.6 dB at 90 Hz and
  then rises to about 81 to 82 dB at 50 to 60 Hz, 5 to 6 dB over its 75.8 dB at 150 Hz (23°, 1/6-octave smoothed,
  `plots/crossover_23_passive.png`). That comes from the bass network's 54 Hz resonance. In a room it is a bass hump.
  Nothing in the fit limits it.
- **The 2 to 5 kHz window.** Presence lift is 4.04 dB passive, 4.63 dB semi-active and 3.90 dB active. All three are
  well above the 0.5 dB the score calls a pass. The lift is the coax's, and it was already there in round 3c, but
  semi-active adds 0.7 dB of it. Window spread is 3.73 dB against 3.40 dB. Above 2851 Hz the scores are low
  confidence and are not ranked.
- **Rear rejection looks better than it is.** The whole_rear score is the minimum over 150 to 300 Hz, so it is set
  near the 150 Hz edge, where passive is a little better. Across the band the picture is different: at about 240 Hz
  the passive rejects about 14 dB, and semi-active and active about 18.5 to 19 dB (`plots/whole_rear_compare.png`).
  The passive's slow filters give less rear rejection in the middle of the band.
- **Sensitivity.** 77.5 dB is 3.3 dB below the woofers alone. It needs about 3.3 dB more amplifier drive than the
  woofers alone for the same level. The active figure is not comparable, as each way has its own amplifier gain.
- **Size and cost.** A 15 mH cored coil rated for 9.6 A, a 20 mH air-core coil, a 120 µF film capacitor, and an
  R_ls that dissipates 86 W at full voltage are not small or cheap. Parts cost was not estimated.

## Notes

- **Models.** The part values come from the model of unmeasured drivers. Real driver impedance will move the
  network. Treat it as a first draft to be trimmed after measurement. The woofers' voice-coil inductance is modelled
  as a pure Le, which overstates their impedance above about 2 kHz.
- **Above 2851 Hz** the data are low confidence. The fit still uses them, at a quarter weight, from 2851 Hz to
  10 kHz, so both coax networks (LF and HF) are shaped partly by low-confidence data. The scores above the limit are
  not ranked.
- **Polarity.** In passive mode the bass and coax HF are wired inverted (driver + terminal to the network's ground
  side), in the wiring between each network and its driver behind K3 and K5. The schematics mark + and − on every
  driver and say "(inverted)". In the active modes polarity is set in the DSP, relative to each driver's own +
  terminal.

## Open

- Whether the owner wants passive at all: it needs one amplifier and no DSP, at lower rear rejection through the cardioid
  band, about 3 dB less sensitivity, a bass hump below 90 Hz, and the parts above. Semi-active needs two channels and
  a DSP on the coax only.
- **Excursion in passive mode.** The resonant coax high-pass gives the coax LF its most voltage at 99 Hz (−5.2 dB
  against the input, more than the −7.5 dB it gets at 1 kHz). The bass network gives the woofers +4.7 dB at 54 Hz.
  Judging either needs the drivers' Xmax and an excursion model.
- The stacked 9.2 dB semi-active coax boost at 139 Hz (with one PEQ at its +12 dB limit) needs a check of the coax
  LF's excursion and amplifier headroom.
- Buy the resistors with margin, or check them against the full-voltage figures in `network.txt`. Check that the
  cored L_b1 is rated for 9.6 A.
- Everything above depends on unmeasured drivers. Rerun after the coax and woofers are measured.
- The VituixCAD netlist (`vituixcad_netlist.txt`) has not been imported or checked by the owner.
- The owner's report page gets its Passive option from the controller after this round.
