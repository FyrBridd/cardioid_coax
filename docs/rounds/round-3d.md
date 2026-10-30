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

**Revised again before merging (owner's choice).** That rerun's bass network peaked at +4.7 dB at 54 Hz, below the
100 Hz the fit looks at, and gave a 5 to 6 dB bass hump at 50 to 60 Hz. The fit now has two guards against that
(item 8 below). The hump is gone, but the guards cost other things on these woofers, and one rule is now broken
(item 8, "What it costs"). Every number below is from the run with the guards.

## What was run

- **Code:** branch `FyrBridd/passive-xo`, commit `5565b51`. The BEM fields came from the cache, so there were no new
  solves. Outputs are not committed.
- **Command:** `uv run speaker passive --params params/round-3.toml --out out/passive/round-3-design --trusted-limit 2851`.
  It took 121 s and exited 0. The log is `out/passive/run.log`. `--trusted-limit` is now required: pass the trusted
  limit of the directivity run.
- **Outputs** (`out/passive/round-3-design/`): `report.html`, `network.txt`, `parts.csv`, `semi_active.txt`,
  `scores.json`, `schematic_bass.svg`, `schematic_coax.svg`, `links.svg`, `vituixcad_netlist.txt` and `plots/`
  (now with `electrical.png`).

## What changed from the plan while building it

1. **The level is tied to the woofers'.** This is a soft target, not a hard rule. The fit's level is the weighted
   mean of the 23° response (1/3-octave smoothed) over 100 Hz to 10 kHz: weight 1 from 100 Hz to the trusted limit
   (2851 Hz), and 0.25 from there to 10 kHz. The fit is penalised when that level falls more than 1 dB under the
   woofers' own level (80-150 Hz). **In this run it sits 1.0 dB under**, inside the 1 dB it may drop for free (the
   run before the bass-peak guards sat 3.5 dB under). The sensitivity is measured differently (23°, 300-1000 Hz),
   and is 0.5 dB under what the woofers give alone. A first preview run (not in this round's outputs) with the level
   left free threw away about 7 dB through a 51 Ω series resistor.
2. **The woofers must stay quiet under the coax LF.** From twice the bass crossover (340 Hz) up to 600 Hz they
   should sit at least 10 dB under the coax LF. This is a penalty in the fit, not a guarantee. Without it, the passive
   woofers ran about 8 dB under the coax out to about 800 Hz (an earlier preview run). That filled the coax LF's
   low-end droop (the cardioid chamber leaves the coax LF about 9 dB down at 150 Hz) and let non-cardioid sound into
   the band the cardioid section must own. **The delivered network nearly meets it, and misses at the bottom of the
   band.** Recomputed from this run (23°, 1/3-octave smoothed, as the fit does): the woofers sit 9.1 dB under the
   coax LF at 356 Hz, so 0.9 dB short, 1 of 11 grid points; 10.0 dB at 378 Hz. From 400 Hz up they are 10.9 dB under
   or more, reaching 15.6 dB at 600 Hz.
3. **A coax LF contour (shelf) network is offered.** It is R_lc with L_lc across it, before the coax LF pad. The run
   before the bass-peak guards used it (5.6 Ω with 0.7 mH, flat23 1.91 to 1.75 dB). **This run leaves it out:** it
   made flat23 worse (2.39 to 2.82 dB). The Zobel across the woofers is used instead (item 8).
4. **Semi-active: the coax DSP channel gets EQ.** The design said gain only. The EQ is fitted before the DSP
   high-pass, as the active fit does, and only inside 150 Hz to 2 times the LF-HF point, cut at the trusted limit
   (150 Hz to 2851 Hz).
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
   4 times past that point on: here the bass from 4 × 170 = 680 Hz up, the coax LF from 6.4 kHz up and below 42 Hz,
   and the coax HF below 400 Hz. The fit penalises every dB over that at 20 per dB (`STOP_WEIGHT`). On this design
   (before item 8) weights 3 and 10 left the bass 0.7 and 1.8 dB over the line near 4 times its crossover point.
   Weights 15 to 100 held it within 0.3 dB. 20 is one step above the lowest that held. `plots/electrical.png` shows
   each branch with the limits.

   **What the rule cost when it was added.** Flatness at 23° went from 1.43 dB (first network) to 1.75 dB. The
   bass-to-coax point moved from 150 Hz to 170 Hz. The level dropped by about 2 dB. It did not force the series
   fallback. **In this run the bass misses it** (item 8): the coax LF and HF stopbands hold (−21.4 dB below 42 Hz,
   −29.1 dB from 6.4 kHz, −23.0 dB below 400 Hz), but the bass is only 5.7 dB under its 170 Hz level at 686 Hz,
   6.3 dB short of the 12 dB.
8. **New, the owner's follow-up: no bass peak below the fitted band.** The run with item 7 gave the woofers
   **+4.7 dB at 54 Hz**: the 15 mH coil resonating with the woofers, which look like a capacitor above their sealed-box
   resonance (about 40 Hz). The sum rose 5 to 6 dB at 50 to 60 Hz. The fit starts at 100 Hz, so nothing limited it.
   It was not fixed by fitting lower: a sealed box rolls off there by design (F3 about 43 Hz, plus room gain), and
   asking it to be flat to 40 Hz invites more resonance. Two guards were added instead:
   - **Electrical gain cap.** No branch may give its driver more than **+1 dB** of voltage (the network on the
     modelled driver loads, 20 Hz-20 kHz). A passive network only gives more than its input by resonating, so this
     holds whatever the acoustic data say.
   - **No acoustic peak below 100 Hz.** Below the fitted band (from 50 Hz, the bottom of the data, up to 100 Hz) the
     23° sum, 1/3-octave smoothed as the fit sees it, may sit at most **1 dB above the fit's level**. A roll-off
     costs nothing.

   Both are penalties, 1000 per dB at each frequency point (`PEAK_WEIGHT`). I tried weights from 20 to 10000 on
   the whole search. **No weight met the +1 dB cap.** In no run did the bass get under +1.1 dB near 60 Hz. My
   reading of the runs: with the 3.5 Ω floor and a bass low-pass at 150 to 190 Hz, a coil small enough not to
   resonate with the woofers leaves too little impedance to hold the floor once a real shunt capacitor is added.
   Weights 20 to 300 left it at +2.0 to +3.6 dB and
   made flat23 3.3 to 7.3 dB. From 500 up it came to +1.1 to +1.7 dB, and the sum stayed within 1 dB of its level
   below 100 Hz. 1000 gave the best flatness, so I kept it. **The result moves a lot with the weight:** 500 gave flat23
   6.06 dB and 2000 gave 1.99 dB, with different crossover points and options.

   **The result.** The woofers now get at most **+1.4 dB, at 61 Hz** (0.4 dB over the cap). The sum below 100 Hz
   rises at most **0.4 dB above the fit's level, at 67 Hz** (limit 1 dB): no hump. Both figures are in `network.txt`
   and the report's summary.

   **What it costs.** To avoid the resonance the fit took a small bass coil (5.5 mH, from 15 mH) and nearly
   dropped the bass capacitor (C_b1 0.56 µF, from 100 µF), and damped the woofers with the Zobel (20 Ω with 390 µF)
   instead. That keeps the floor, but the bass low-pass is now gentle: the woofers get −13.5 dB at 680 Hz, −14 dB at
   1 kHz and −15 to −16 dB from 2 to 3 kHz. **So the bass misses the stopband rule of item 7 by 6.3 dB.** Every
   weight I checked it at (50, 150, 1000, 2000, 5000 and 10000) missed it, by 6 to 18 dB. Against the run before
   the guards: flat23 1.75 to **1.42 dB** (better), sensitivity 77.5 to **80.3 dB** (better), minimum impedance
   3.52 to **3.50 Ω** (still on the floor), but whole-system rear rejection 11.62 to **7.93 dB**, bass DI smoothness
   0.69 to **1.32 dB**, presence lift 4.04 to **4.32 dB** and window spread 3.73 to **3.99 dB** (all worse). The
   crossover points moved from 170 / 1200 Hz to **170 / 1600 Hz**. It does not force the series fallback.

## The passive network

**Crossover points: bass to coax 170 Hz, coax LF to HF 1600 Hz.** The search took 30 grid pairs and kept these.
It uses the woofers in parallel, the bass inverted and the coax HF inverted. One amplifier, one pair of terminals
(IN), and links K1 to K5 all closed. Every value is a standard value (capacitors E12, resistors E24).

| Part | Value | Type | DCR | Peak current at full voltage | Resistor rating |
|---|---|---|---|---|---|
| L_b1 | 5.5 mH | cored inductor | 0.33 Ω | 10.1 A (saturation current at least this) | |
| C_b1 | 0.56 µF | film capacitor | | | |
| R_bz | 20 Ω | resistor (Zobel) | | | 5 W |
| C_bz | 390 µF | bipolar electrolytic capacitor (Zobel) | | | |
| C_h1 | 47 µF | film capacitor | | | |
| L_h1 | 20 mH | air-core inductor | 2.24 Ω | 1.9 A | |
| L_l1 | 1.5 mH | air-core inductor | 0.61 Ω | 6.5 A | |
| C_l1 | 10 µF | film capacitor | | | |
| R_ls | 20 Ω | resistor | | | 10 W |
| R_lp | 22 Ω | resistor | | | 5 W |
| C_t1 | 4.7 µF | film capacitor | | | |
| L_t1 | 0.5 mH | air-core inductor | 0.35 Ω | 1.1 A | |
| C_t2 | 390 µF | film capacitor | | | |
| R_ts | 2.2 Ω | resistor | | | 15 W |
| R_tp | 1.2 Ω | resistor | | | 10 W |

The peak current is the peak of a sine with the amplifier's full 25 V rms on the input, the largest from 20 Hz to
20 kHz in either mode. A cored coil must not saturate below it. Parts per pair are in `parts.csv` (column
`peak_a`).

**What each part does.**
- **Bass (L_b1, C_b1, and the Zobel R_bz, C_bz).** A series coil, and across the woofers a 20 Ω resistor in series
  with 390 µF. C_b1 is so small (0.56 µF, about 1.7 kΩ at 170 Hz) that it does almost nothing; the Zobel arm is
  the real shunt, and its 20 Ω damps the coil against the woofers. Its transfer (`plots/electrical.png`) rises to
  **+1.4 dB at 61 Hz** (it was +4.7 dB at 54 Hz before the guards), is −2.6 dB at 100 Hz and −7.8 dB at 170 Hz, then
  falls gently: −13.5 dB at 680 Hz, −14 dB at 1 kHz, −15.1 dB at 2 kHz, −16.3 dB at 3 kHz and −19.1 dB at 5 kHz. That
  is a 1st-order slope into a shelf, not a 2nd-order low-pass, and it misses the stopband rule by 6.3 dB at 686 Hz
  (item 8).
- **Coax high-pass (C_h1 from IN to node cx, L_h1 from cx to ground).** It feeds both coax networks, the LF and the
  HF. The coax LF gets −15.7 dB at 99 Hz and −11.0 dB at 150 Hz.
- **Coax LF.** A low-pass (L_l1, C_l1) and a pad (R_ls, R_lp), no contour. It is resonant near its crossover: the
  coax LF gets its most voltage, −5.0 dB, at 1.41 kHz.
- **Coax HF.** A high-pass (C_t1, L_t1, C_t2) and a pad (R_ts, R_tp). C_t2 (390 µF) is so large that it is nearly a
  wire, so the high-pass is close to 2nd order. A 390 µF film capacitor is a big, costly part; the model has not been
  run without it.

**What the search kept and left out:**
- Kept the Zobel across the woofers (flat23 2.39 to 1.42 dB).
- Left out the coax LF contour (flat23 2.39 to 2.82 dB), a 3rd-order coax LF low-pass (1.42 to 2.86 dB) and the RLC
  shunt across the coax HF (1.42 to 1.99 dB): each made flat23 worse.
- Left out a 3rd-order bass low-pass: its fit took the impedance to 0.15 Ω, far under the floor.
- Left out the coax LF notch at 2700 Hz: it made presence lift 0.40 dB worse (0.25 dB allowed).

**Impedance and level:**
- Minimum impedance is **3.50 Ω at 73 Hz**, on the 3.5 Ω floor (+0.00 Ω). That is also the minimum in the bass band
  (20-500 Hz). `plots/impedance.png` shows the 3.5 Ω line. The floor is the SA550's minimum.
- **Semi-active impedance.** The coax channel's minimum is 3.65 Ω at 20 kHz, 0.15 Ω over the floor. The bass channel
  (bare woofers in parallel) dips to **3.26 Ω at 149 Hz, below the 3.5 Ω the SA550 needs**. The fit and the tolerance
  check do not cover that channel, because the amplifier for it is not chosen. Pick one that is happy at about 3 Ω.
- Sensitivity is **80.3 dB at 2 m for 2.83 V** (23°, 300-1000 Hz). That is 0.5 dB below what the woofers give alone
  (77.5 dB, 3.3 dB below, before the guards).
- **Tolerance check passed.** The worst of 200 draws was 3.47 Ω (limit 3.3 Ω). Flatness (flat23) was 1.42 dB nominal
  and 1.67 dB at the 95th percentile (limit +0.5 dB). Whole_rear ran 7.57 to 8.21 dB (5th to 95th percentile).
- No series-woofer fallback was needed: a parallel network cleared the floor.

**Power.** Resistor ratings follow the rule in `network.txt`: the first of 5, 10, 15, 20, 30, 50 W at or above the
dissipation with one third of the amplifier's voltage (roughly what music averages). Three are close to their
ratings: R_bz at 4.31 W against 5 W (passive mode only), and, in semi-active mode (the same coax boards, but the
amplifier channel carries only the coax), R_ts at 11.61 W against 15 W and R_ls at 7.63 W against 10 W (passive
mode: 11.18 W and 7.30 W). R_tp reaches 5.42 W against 10 W. At full voltage, in the larger of the two modes, R_ts
would dissipate 105 W, R_ls 69 W, R_tp 49 W and R_bz 39 W. I would buy higher-rated parts than the listed ones, at
least for R_bz and R_ts. These are sine figures, so they are conservative for music.

**Coils and capacitors.** L_b1 is 5.5 mH on a core and must carry 10.1 A peak without saturating. L_h1 (20 mH air
core, 2.24 Ω) is large. C_bz (390 µF electrolytic) and C_t2 (390 µF film) are large too.

### Warnings in network.txt

- **C_bz is at the 400 µF limit.** The fit wanted a larger Zobel capacitor. The response near this part is a
  compromise.
- **L_h1 is at the 20 mH limit.** The fit wanted an even larger coil for the coax high-pass, but 20 mH air-core is
  already 2.24 Ω of DCR and would be large and costly. The response near this part is a compromise.
- `network.txt` also marks two lines **OVER**: the bass stopband (−5.7 dB at 686 Hz, limit −12 dB) and the largest
  branch gain (+1.4 dB at 61 Hz, limit +1 dB). See item 8.

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
| bass | +0.10 dB | 1.30 ms | inverted | 80 Hz +1.18 dB Q1.65; 214 Hz −1.77 dB Q1.22 |
| coax | 0 dB | 0 | normal | 150 Hz +8.45 dB Q1.27; 598 Hz −2.46 dB Q2.25; 1322 Hz −7.09 dB Q2.00; 1449 Hz +6.24 dB Q1.20 |

- **No PEQ is at the fitter's limit now.** The coax boards changed with the bass-peak guards (no contour, a lighter
  coax LF pad), and in semi-active mode they pass the coax LF −10.8 dB at 150 Hz, about 4 dB more than the boards
  before the guards. So the DSP has less of the cardioid chamber's low-end droop to fill: +8.45 dB at 150 Hz instead
  of +12 dB at its limit.
- **A cut and a boost close together.** The 1322 Hz cut and the 1449 Hz boost partly cancel. Together they make
  about −1.2 dB at 1.3 kHz and +2 dB at 1.8 to 2 kHz, around the coax LF to HF crossover.
- **Stacked boost.** The four coax PEQs together (my sum of their curves) peak at **8.4 dB at 150 Hz**. That needs
  the usual care: amplifier clipping and the coax LF's excursion near 150 Hz should be checked (see Open). Normalise in
  the DSP so the largest boost does not clip, as the file says.

## Three modes, trusted scores

Trusted band (150 Hz to 2851 Hz; bass DI from 120 Hz). Lower is better except the rear and cardioid rows. The first
two columns are earlier passive networks, for comparison: the first network, and the one with the stopband rule but
before the bass-peak guards. Semi-active and active are from this run.

| Score (dB) | Passive, first network | Passive, before the guards | Passive | Semi-active | Active |
|---|---|---|---|---|---|
| 23° flatness | 1.43 | 1.75 | 1.42 | **0.81** | 1.16 |
| Window spread | 3.60 | 3.73 | 3.99 | 3.97 | **3.40** |
| Presence lift, 2-2.85 kHz | 4.22 | 4.04 | 4.32 | 4.23 | **3.90** |
| DI smoothness, bass, 120-500 Hz | 0.94 | 0.69 | 1.32 | 1.53 | 1.66 |
| DI smoothness, LF-HF, 500-2851 Hz | 0.77 | 1.11 | 1.20 | 1.12 | **1.08** |
| Whole-system rear rejection, 150-300 Hz (higher is better) | 10.44 | 11.62 | 7.93 | 9.87 | **9.88** |
| Cardioid hard / aim (coax section) | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 |

Bold marks the best of this run's three modes. The coax section's cardioid scores are the same because the coax is
the same box.

## Where passive falls short, and why

- **Flatness.** Passive is 0.26 dB behind active at 23°. Semi-active is now 0.35 dB ahead of active: its DSP EQ
  works on the new coax boards without hitting a limit. The owner can judge all three from
  `plots/crossover_23_*.png`.
- **Below the fitted band: no hump now.** Below 100 Hz the passive sum (23°, 1/6-octave smoothed,
  `plots/crossover_23_passive.png`) is 78.2 dB at 50 Hz, 79.7 dB at 60 Hz and 80.1 dB at 67 to 70 Hz, then dips to
  76.9 dB at 120 Hz and comes back to 78.8 dB at 150 Hz. Against the fit's level (1/3-octave, as the fit sees it) the
  highest point below 100 Hz is 0.4 dB above, at 67 Hz. Before the guards it was 5 to 6 dB over its 150 Hz level at
  50 to 60 Hz.
- **The bass reaches into the coax band.** The bass low-pass is now gentle (item 8): the woofers get −14 to −16 dB of
  the input from 1 to 3 kHz, against −34 dB at 1 kHz before. The model shows little sound from them there, but the
  bass field is not modelled above the trusted limit (2851 Hz), and real woofers have cone break-up in that region.
  With the owner sensitive to 2 to 5 kHz, this needs a measurement before building passive.
- **The 2 to 5 kHz window.** Presence lift is 4.32 dB passive, 4.23 dB semi-active and 3.90 dB active. All three are
  well above the 0.5 dB the score calls a pass. The lift is the coax's, and it was already there in round 3c;
  passive adds 0.4 dB of it and semi-active 0.3 dB. Window spread is 3.99 and 3.97 dB against 3.40 dB. Above 2851 Hz
  the scores are low confidence and are not ranked.
- **Rear rejection is lower.** The whole_rear score is the minimum over 150 to 300 Hz. Passive is 7.93 dB against
  9.87 and 9.88 dB, and across the band it is lower too: at about 240 Hz the passive rejects about 13.4 dB, and
  semi-active and active about 18.5 dB (`plots/whole_rear_compare.png`). The gentle bass low-pass lets the woofers,
  which radiate both ways, further into the cardioid band.
- **Sensitivity.** 80.3 dB is 0.5 dB below the woofers alone. The active figure is not comparable, as each way has
  its own amplifier gain.
- **Size and cost.** A 5.5 mH cored coil rated for 10.1 A, a 20 mH air-core coil, a 390 µF film capacitor, a
  390 µF electrolytic, and an R_ts that dissipates 105 W at full voltage are not small or cheap. Parts cost was not
  estimated.

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

- Whether the owner wants passive at all: it needs one amplifier and no DSP, at lower rear rejection through the
  cardioid band, a bass low-pass that misses its stopband rule and reaches into the coax band, and the parts above.
  Semi-active needs two channels and a DSP on the coax only, and is now the flattest of the three.
- **The rules conflict on these woofers.** The +1 dB gain cap, the 3.5 Ω floor and the bass stopband rule cannot all
  hold with the topologies offered: no weight tried met the cap, and every weight that came near it broke the
  stopband rule (item 8). Ways out, none modelled here: accept the gentle bass slope (this run), accept the 54 Hz
  peak (the run before), or add a network that flattens the woofers' impedance peak (a series RLC across them, tuned
  to their resonance), which would let a normal 2nd-order low-pass work. The owner chooses.
- **The result is sensitive.** Small changes to the guard weight changed the crossover points, the options and
  flat23 a lot (item 8). Measured driver data may move it as much.
- **Excursion in passive mode.** The woofers get +1.4 dB at 61 Hz. The coax LF gets its most voltage, −5.0 dB, at
  1.41 kHz, and −15.7 dB at 99 Hz. Judging excursion needs the drivers' Xmax and an excursion model.
- The stacked 8.4 dB semi-active coax boost at 150 Hz needs a check of the coax LF's excursion and amplifier
  headroom.
- Buy the resistors with margin, or check them against the full-voltage figures in `network.txt`. Check that the
  cored L_b1 is rated for 10.1 A.
- Everything above depends on unmeasured drivers. Rerun after the coax and woofers are measured.
- The VituixCAD netlist (`vituixcad_netlist.txt`) has not been imported or checked by the owner.
- The owner's report page gets its Passive option from the controller after this round.
