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
(item 8 below). The hump went, but on these woofers the guards alone broke the bass stopband rule (item 8).

**Revised a third time (owner chose the compensator).** The search now offers an impedance-peak compensator across
the woofers (item 9). With it every rule holds: the floor, the stopband rule, the gain cap and the no-hump guard.
Every number below is from the run with the compensator.

## What was run

- **Code:** branch `FyrBridd/passive-xo`, commit `d1a4a7b`. The BEM fields came from the cache, so there were no new
  solves. Outputs are not committed.
- **Command:** `uv run speaker passive --params params/round-3.toml --out out/passive/round-3-design --trusted-limit 2851`.
  It took about 2 minutes and exited 0. The log is `out/passive/run.log`. `--trusted-limit` is now required: pass the
  trusted limit of the directivity run.
- **Outputs** (`out/passive/round-3-design/`): `report.html`, `network.txt`, `parts.csv`, `semi_active.txt`,
  `scores.json`, `schematic_bass.svg`, `schematic_coax.svg`, `links.svg`, `vituixcad_netlist.txt` and `plots/`
  (now with `electrical.png`).

## What changed from the plan while building it

1. **The level is tied to the woofers'.** This is a soft target, not a hard rule. The fit's level is the weighted
   mean of the 23° response (1/3-octave smoothed) over 100 Hz to 10 kHz: weight 1 from 100 Hz to the trusted limit
   (2851 Hz), and 0.25 from there to 10 kHz. The fit is penalised when that level falls more than 1 dB under the
   woofers' own level (80-150 Hz). **In this run it sits 1.1 dB under**, 0.1 dB past the 1 dB it may drop for free
   (the run with the guards alone sat 1.0 dB under, the run before them 3.5 dB). The sensitivity is measured
   differently (23°, 300-1000 Hz), and is 1.4 dB under what the woofers give alone. A first preview run (not in this
   round's outputs) with the level left free threw away about 7 dB through a 51 Ω series resistor.
2. **The woofers must stay quiet under the coax LF.** From twice the bass crossover (340 Hz) up to 600 Hz they
   should sit at least 10 dB under the coax LF. This is a penalty in the fit, not a guarantee. Without it, the passive
   woofers ran about 8 dB under the coax out to about 800 Hz (an earlier preview run). That filled the coax LF's
   low-end droop (the cardioid chamber leaves the coax LF about 9 dB down at 150 Hz) and let non-cardioid sound into
   the band the cardioid section must own. **The delivered network nearly meets it, and misses at the bottom of the
   band.** Recomputed from this run (23°, 1/3-octave smoothed, as the fit does): the woofers sit 9.7 dB under the
   coax LF at 356 Hz, so 0.3 dB short, 1 of 11 grid points; 12.0 dB at 378 Hz. From 400 Hz up they are 14.9 dB under
   or more, reaching 27.1 dB at 600 Hz.
3. **A coax LF contour (shelf) network is offered.** It is R_lc with L_lc across it, before the coax LF pad. The run
   before the bass-peak guards used it (5.6 Ω with 0.7 mH, flat23 1.91 to 1.75 dB). **This run leaves it out:** it
   made flat23 a little worse (1.62 to 1.67 dB).
4. **Semi-active: the coax DSP channel gets EQ.** The design said gain only. The EQ is fitted before the DSP
   high-pass, as the active fit does, and only inside 150 Hz to 2 times the LF-HF point, cut at the trusted limit
   (150 Hz to 2851 Hz).
5. **Terminal links K1 to K5.** Each network is isolated at its boundary, so the same boards serve all three modes
   (`links.svg`, and the table at the end of `network.txt`). **K3, K4 and K5 are two-pole:** each breaks both the hot
   and the return between a network and its driver. In the active modes every terminal pair then reaches its driver
   alone. Without that, the drivers' returns stayed joined at the network ground, and in passive mode the inverted
   woofers have their + on that ground: wiring the BASS terminals "straight to the woofers" could have shorted an
   amplifier channel. The compensator (item 9) sits on the network side of K3, like the Zobel, so it drops out with
   the bass network in the active modes.
6. **Impedance and tolerance rules.** The floor is 3.5 Ω. The tolerance check draws every inductor and capacitor at
   plus or minus 5 % (200 draws). If it fails, the fit is redone with the floor 0.1 Ω higher (up to 3.8 Ω). **In this
   run it failed at 3.5 Ω** (worst draw 3.26 Ω, limit 3.3 Ω), **and passed at 3.6 Ω.** The series-woofer fallback is
   tried only if no parallel network clears 3.5 Ω.
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
   fallback. **In this run every stopband holds:** the bass is 21.0 dB under its 170 Hz level at 686 Hz (it missed by
   6.3 dB with the guards alone), the coax LF −14.1 dB below 42 Hz and −28.7 dB from 6.4 kHz, the coax HF −22.8 dB
   below 400 Hz.
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
   the whole search, before the compensator existed. **No weight met the +1 dB cap.** In no run did the bass get
   under +1.1 dB near 60 Hz. With the 3.5 Ω floor and a bass low-pass at 150 to 190 Hz, a coil small enough not to
   resonate with the woofers left too little impedance to hold the floor once a real shunt capacitor was added.
   Weights 20 to 300 left it at +2.0 to +3.6 dB and made flat23 3.3 to 7.3 dB. From 500 up it came to +1.1 to
   +1.7 dB. 1000 gave the best flatness, so I kept it. **The result moved a lot with the weight:** 500 gave flat23
   6.06 dB and 2000 gave 1.99 dB, with different crossover points and options.

   **What the guards alone cost** (the run before item 9). The woofers got +1.4 dB at 61 Hz, 0.4 dB over the cap. To
   hold that, the fit took a small bass coil (5.5 mH) and nearly dropped the bass capacitor (0.56 µF), with a Zobel
   (20 Ω with 390 µF) as the real shunt. The bass low-pass became a gentle slope into a shelf (−13.5 dB at 680 Hz,
   −15 to −16 dB from 2 to 3 kHz) and **missed the stopband rule of item 7 by 6.3 dB.** Whole-system rear rejection
   fell to 7.93 dB. The owner chose to fix that with the compensator (item 9).
9. **New, the owner's choice: an impedance-peak compensator across the woofers.** The woofer pair in its sealed box
   has an impedance peak of about **29 Ω at 40 Hz** (the box resonance), and falls to 3.26 Ω at 147 Hz. Above the
   peak the woofers look like a capacitor, and that is what the bass coil resonated with (item 8). The compensator is
   the standard fix: a resistor, a coil and a capacitor in series (R_bc, L_bc, C_bc), wired across the woofers on the
   network side of K3. Near the box resonance it takes current the woofers do not, so the network sees a much flatter
   load and a normal low-pass can work into it. Its textbook values come from the woofers' Re, fc, Qec and Qmc:
   3.6 Ω, 9.2 mH and 1600 µF here. With those, the pair's 20-200 Hz impedance would run from 2.7 to 3.3 Ω instead of
   3.3 to 29 Ω. That is under the 3.5 Ω floor, so the fit detunes it: **2.7 Ω, 23 mH and 1200 µF**, tuned to about
   30 Hz. With these the woofer load runs from 2.6 Ω at 30 Hz to 6.2 Ω at 55 Hz (31 Ω without), and the network's
   series coils lift the input to the floor. It is the first option the search tries, and it was kept: **flat23 2.39
   to 1.62 dB.** After it, the search kept a 3rd-order bass low-pass (1.62 to 1.45 dB), and left the Zobel out
   (1.53 dB, less than the 0.1 dB gain it needs).

   **The result.** The woofers now get at most **−0.3 dB, at 64 Hz** (limit +1 dB). The sum below 100 Hz stays
   **1.2 dB under the fit's level** at its highest (67 Hz; limit +1 dB above it): no hump. The bass meets the
   stopband rule with 9 dB to spare (item 7). Whole-system rear rejection is back to **9.70 dB**.

   **What it costs.** Three more parts, two of them large: L_bc is a 23 mH cored coil that must carry 7.7 A peak
   without saturating, and C_bc a 1200 µF bipolar electrolytic. R_bc (2.7 Ω) dissipates 8.9 W with a third of the
   amplifier's voltage (10 W part) and 80 W at full voltage. The bass low-pass also grew: L_b1 8 mH (10.0 A) and L_b2
   0.5 mH (10.1 A), both cored, and C_b1 180 µF. Against the run with the guards alone: flat23 1.42 to **1.68 dB**,
   sensitivity 80.3 to **79.4 dB**, presence lift 4.32 to **4.72 dB** and window spread 3.99 to **4.20 dB** (all
   worse); bass DI smoothness 1.32 to **0.74 dB** and whole_rear 7.93 to **9.70 dB** (better). The coax boards
   changed too, which moves semi-active: its flat23 went 0.81 to 1.16 dB and presence lift 4.23 to 4.55 dB. The
   crossover points stay at 170 / 1600 Hz. No series fallback.

## The passive network

**Crossover points: bass to coax 170 Hz, coax LF to HF 1600 Hz.** The search took 30 grid pairs and kept these.
It uses the woofers in parallel, the bass inverted and the coax HF inverted. One amplifier, one pair of terminals
(IN), and links K1 to K5 all closed. Every value is a standard value (capacitors E12, resistors E24). The fit's
impedance floor is 3.6 Ω (item 6).

| Part | Value | Type | DCR | Peak current at full voltage | Resistor rating |
|---|---|---|---|---|---|
| L_b1 | 8 mH | cored inductor | 0.40 Ω | 10.0 A (saturation current at least this) | |
| C_b1 | 180 µF | bipolar electrolytic capacitor | | | |
| L_b2 | 0.5 mH | cored inductor | 0.10 Ω | 10.1 A (saturation current at least this) | |
| R_bc | 2.7 Ω | resistor (compensator) | | | 10 W |
| L_bc | 23 mH | cored inductor (compensator) | 0.67 Ω | 7.7 A (saturation current at least this) | |
| C_bc | 1200 µF | bipolar electrolytic capacitor (compensator) | | | |
| C_h1 | 120 µF | film capacitor | | | |
| L_h1 | 20 mH | air-core inductor | 2.24 Ω | 4.7 A | |
| L_l1 | 1.5 mH | air-core inductor | 0.61 Ω | 6.1 A | |
| C_l1 | 8.2 µF | film capacitor | | | |
| R_ls | 24 Ω | resistor | | | 10 W |
| R_lp | 18 Ω | resistor | | | 5 W |
| C_t1 | 1 µF | film capacitor | | | |
| L_t1 | 0.7 mH | air-core inductor | 0.42 Ω | 0.34 A | |
| C_t2 | 390 µF | film capacitor | | | |
| R_ts | 0.68 Ω | resistor | | | 5 W |
| R_tp | 12 Ω | resistor | | | 5 W |

The peak current is the peak of a sine with the amplifier's full 25 V rms on the input, the largest from 20 Hz to
20 kHz in either mode. A cored coil must not saturate below it. Parts per pair are in `parts.csv` (column
`peak_a`).

**What each part does.**
- **Bass (L_b1, C_b1, L_b2, and the compensator R_bc, L_bc, C_bc).** A 3rd-order low-pass (coil, capacitor to
  ground, coil) into the woofers, with the compensator across them (item 9). Its transfer (`plots/electrical.png`)
  is −3.3 dB at 35 Hz, at most **−0.3 dB at 64 Hz**, −4.0 dB at 100 Hz and −8.2 dB at 170 Hz. It holds about −9 to
  −10 dB from 200 to 300 Hz, then falls steeply: −29 dB at 680 Hz, −37 dB at 1 kHz, −50 dB at 2 kHz. That meets the
  stopband rule (−21.0 dB re 170 Hz at 686 Hz, limit −12 dB).
- **Coax high-pass (C_h1 from IN to node cx, L_h1 from cx to ground).** It feeds both coax networks, the LF and the
  HF. It is now resonant near 100 Hz: the coax LF gets **−6.5 dB at 100 Hz** (it was −15.7 dB at 99 Hz with the
  guards alone) and −9.6 dB at 150 Hz.
- **Coax LF.** A low-pass (L_l1, C_l1) and a pad (R_ls, R_lp), no contour. It is resonant near its crossover: the
  coax LF gets its most voltage, −6.0 dB, at 1.41 kHz.
- **Coax HF.** C_t1 (1 µF) with L_t1 (0.7 mH) is a high-pass tuned near 6 kHz, which here does most of the padding
  as well: the HF gets −24 dB at 2 kHz, −20 dB at 3 kHz, −16 dB at 5 kHz and −4 dB at 20 kHz, rising with frequency.
  The pad (R_ts 0.68 Ω, R_tp 12 Ω) does little. C_t2 (390 µF) is so large that it is nearly a wire; a 390 µF film
  capacitor is a big, costly part, and the model has not been run without it.

**What the search kept and left out** (each option tried in turn, in this order):
- Kept the impedance-peak compensator across the woofers (flat23 2.39 to 1.62 dB).
- Left out the coax LF contour (1.62 to 1.67 dB): no better.
- Left out the Zobel across the woofers (1.62 to 1.53 dB): less than the 0.1 dB gain needed.
- Kept a 3rd-order bass low-pass (1.62 to 1.45 dB).
- Left out a 3rd-order coax LF low-pass (1.45 to 1.62 dB) and the coax LF notch at 2700 Hz (1.45 to 1.47 dB): no
  better.
- Left out the RLC shunt across the coax HF: it made the LF-HF directivity smoothness 1.37 dB worse (0.25 dB
  allowed).

The tolerance refit at the 3.6 Ω floor moved flat23 from 1.45 to 1.68 dB.

**Impedance and level:**
- Minimum impedance is **3.60 Ω at 25 Hz**, on the 3.6 Ω floor the fit used, and 0.10 Ω over the 3.5 Ω the SA550
  needs. That is also the minimum in the bass band (20-500 Hz). `plots/impedance.png` shows the 3.5 Ω line.
- **Semi-active impedance.** The coax channel's minimum is 6.27 Ω at 1436 Hz. The bass channel (bare woofers in
  parallel, no compensator: it is behind K3) dips to **3.26 Ω at 149 Hz, below the 3.5 Ω the SA550 needs**. The fit
  and the tolerance check do not cover that channel, because the amplifier for it is not chosen. Pick one that is
  happy at about 3 Ω.
- Sensitivity is **79.4 dB at 2 m for 2.83 V** (23°, 300-1000 Hz). That is 1.4 dB below what the woofers give alone
  (80.3 dB with the guards alone, 77.5 dB before them).
- **Tolerance check passed at 3.6 Ω.** The worst of 200 draws was 3.41 Ω (limit 3.3 Ω). Flatness (flat23) was
  1.68 dB nominal and 1.91 dB at the 95th percentile (limit +0.5 dB). Whole_rear ran 9.48 to 9.93 dB (5th to 95th
  percentile).
- No series-woofer fallback was needed: a parallel network cleared the floor.

**Power.** Resistor ratings follow the rule in `network.txt`: the first of 5, 10, 15, 20, 30, 50 W at or above the
dissipation with one third of the amplifier's voltage (roughly what music averages). Two are close to their ratings:
R_bc at 8.87 W against 10 W (passive mode only), and R_ls at 7.97 W against 10 W (7.61 W in semi-active mode). R_tp
reaches 2.91 W against 5 W, R_lp 0.98 W and R_ts 0.57 W. At full voltage, in the larger of the two modes, R_bc would
dissipate 80 W, R_ls 72 W, R_tp 26 W, R_lp 9 W and R_ts 5 W. I would buy higher-rated parts than the listed ones,
at least for R_bc and R_ls. These are sine figures, so they are conservative for music.

**Coils and capacitors.** Three cored coils: L_b1 (8 mH, 10.0 A peak), L_b2 (0.5 mH, 10.1 A) and L_bc (23 mH, 7.7 A),
each rated for at least that saturation current. L_h1 (20 mH air core, 2.24 Ω) is large. C_bc (1200 µF) and C_b1
(180 µF) are bipolar electrolytics; C_h1 (120 µF) and C_t2 (390 µF) are large film capacitors.

### Warnings in network.txt

- **L_h1 is at the 20 mH limit.** The fit wanted an even larger coil for the coax high-pass, but 20 mH air-core is
  already 2.24 Ω of DCR and would be large and costly. The response near this part is a compromise.
- No line is marked **OVER**: the stopbands, the largest branch gain (−0.3 dB at 65 Hz, limit +1 dB) and the rise
  below 100 Hz (−1.2 dB at 67 Hz, limit +1 dB) all hold.

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
| bass | −1.79 dB | 1.30 ms | inverted | 80 Hz +1.18 dB Q1.65; 214 Hz −1.77 dB Q1.22 |
| coax | 0 dB | 0 | normal | 150 Hz +8.28 dB Q1.32; 572 Hz −2.07 dB Q3.13; 1352 Hz −2.07 dB Q5.13; 2726 Hz −1.15 dB Q3.27 |

- **No PEQ is at the fitter's limit.** In semi-active mode the coax boards pass the coax LF −9.4 dB at 150 Hz, so the
  DSP fills the cardioid chamber's low-end droop with +8.28 dB at 150 Hz.
- **Stacked boost.** The four coax PEQs together (my sum of their curves) peak at **8.3 dB at 150 Hz**. That needs
  the usual care: amplifier clipping and the coax LF's excursion near 150 Hz should be checked (see Open). Normalise in
  the DSP so the largest boost does not clip, as the file says.

## Three modes, trusted scores

Trusted band (150 Hz to 2851 Hz; bass DI from 120 Hz). Lower is better except the rear and cardioid rows. The first
three columns are earlier passive networks, for comparison: the first network, the one with the stopband rule but
before the bass-peak guards, and the one with the guards but no compensator. Semi-active and active are from this
run.

| Score (dB) | Passive, first network | Passive, before the guards | Passive, guards alone | Passive | Semi-active | Active |
|---|---|---|---|---|---|---|
| 23° flatness | 1.43 | 1.75 | 1.42 | 1.68 | **1.16** | **1.16** |
| Window spread | 3.60 | 3.73 | 3.99 | 4.20 | 4.20 | **3.40** |
| Presence lift, 2-2.85 kHz | 4.22 | 4.04 | 4.32 | 4.72 | 4.55 | **3.90** |
| DI smoothness, bass, 120-500 Hz | 0.94 | 0.69 | 1.32 | **0.74** | 1.53 | 1.66 |
| DI smoothness, LF-HF, 500-2851 Hz | 0.77 | 1.11 | 1.20 | 1.25 | 1.24 | **1.08** |
| Whole-system rear rejection, 150-300 Hz (higher is better) | 10.44 | 11.62 | 7.93 | 9.70 | **9.89** | 9.88 |
| Cardioid hard / aim (coax section) | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 |

Bold marks the best of this run's three modes. The coax section's cardioid scores are the same because the coax is
the same box. Semi-active's scores moved with the coax boards (with the guards alone it had flat23 0.81, window
spread 3.97, presence lift 4.23, DI LF-HF 1.12 and whole_rear 9.87 dB).

## Where passive falls short, and why

- **Flatness.** Passive is 0.52 dB behind active and semi-active at 23° (1.68 against 1.16 dB). The owner can judge
  all three from `plots/crossover_23_*.png`.
- **Below the fitted band: no hump.** Below 100 Hz the passive sum (23°, 1/6-octave smoothed,
  `plots/crossover_23_passive.png`) is 75.7 dB at 50 Hz, 77.7 dB at 60 Hz, 78.5 dB at 67 Hz and 77.6 dB at 100 Hz,
  then rises to 79.4 dB at 150 Hz and 81.0 dB at 250 Hz. Against the fit's level (1/3-octave, as the fit sees it)
  the highest point below 100 Hz is 1.2 dB under, at 67 Hz.
- **The 2 to 5 kHz window.** Presence lift is 4.72 dB passive, 4.55 dB semi-active and 3.90 dB active. All three are
  well above the 0.5 dB the score calls a pass. The lift is the coax's, and it was already there in round 3c; passive
  adds 0.8 dB of it and semi-active 0.65 dB, more than with the guards alone (0.4 and 0.3 dB). The coax HF network
  changed (item 9 and "What each part does"), and that is the likely cause. Window spread is 4.20 dB in both against
  3.40 dB. Above 2851 Hz the scores are low confidence and are not ranked. The bass no longer reaches into this band:
  the woofers get −50 dB at 2 kHz (−15 dB with the guards alone).
- **Rear rejection.** The whole_rear score is the minimum over 150 to 300 Hz. Passive is 9.70 dB against 9.89 and
  9.88 dB, close now, but across the band it is lower: at about 240 Hz the passive rejects about 12.8 dB, and
  semi-active and active about 18.5 dB (`plots/whole_rear_compare.png`). Below 150 Hz the passive rejects more.
- **Sensitivity.** 79.4 dB is 1.4 dB below the woofers alone. The active figure is not comparable, as each way has
  its own amplifier gain.
- **Size and cost.** Three cored coils (8 mH and 0.5 mH at 10 A, 23 mH at 7.7 A), a 20 mH air-core coil, a 1200 µF
  and a 180 µF electrolytic, a 120 µF and a 390 µF film capacitor, and an R_bc that dissipates 80 W at full voltage
  are not small or cheap. Parts cost was not estimated.

## Notes

- **Models.** The part values come from the model of unmeasured drivers. Real driver impedance will move the
  network, and the compensator most of all: it is tuned to the woofers' modelled box resonance. Treat it as a first
  draft to be trimmed after measurement. The woofers' voice-coil inductance is modelled as a pure Le, which
  overstates their impedance above about 2 kHz.
- **Above 2851 Hz** the data are low confidence. The fit still uses them, at a quarter weight, from 2851 Hz to
  10 kHz, so both coax networks (LF and HF) are shaped partly by low-confidence data. The scores above the limit are
  not ranked.
- **Polarity.** In passive mode the bass and coax HF are wired inverted (driver + terminal to the network's ground
  side), in the wiring between each network and its driver behind K3 and K5. The schematics mark + and − on every
  driver and say "(inverted)". In the active modes polarity is set in the DSP, relative to each driver's own +
  terminal.

## Open

- Whether the owner wants passive at all: it needs one amplifier and no DSP, at 0.5 dB less flatness than the other
  two, more 2 to 5 kHz lift, lower rear rejection through the middle of the cardioid band, and the parts above.
  Semi-active needs two channels and a DSP on the coax only, and ties active for flatness.
- **The rule conflict is resolved by the compensator.** With it the floor, the gain cap, the no-hump guard and the
  bass stopband rule all hold (item 9). The cost is three more parts and the coax changes above.
- **The result is sensitive.** Small changes to the guard weight changed the crossover points, the options and
  flat23 a lot (item 8), and the compensator's values depend on the woofers' modelled resonance. Measured driver
  data may move it as much.
- **Excursion in passive mode.** The woofers get at most −0.3 dB, at 64 Hz. The coax LF now gets −6.5 dB at 100 Hz,
  9 dB more than with the guards alone, and its most, −6.0 dB, at 1.41 kHz. Judging excursion needs the drivers'
  Xmax and an excursion model; the coax LF near 100 Hz is the one to check.
- The stacked 8.3 dB semi-active coax boost at 150 Hz needs a check of the coax LF's excursion and amplifier
  headroom.
- Buy the resistors with margin, or check them against the full-voltage figures in `network.txt`. Check that the
  cored L_b1, L_b2 and L_bc are rated for 10.0, 10.1 and 7.7 A.
- Everything above depends on unmeasured drivers. Rerun after the coax and woofers are measured.
- The VituixCAD netlist (`vituixcad_netlist.txt`) has not been imported or checked by the owner.
- The owner's report page gets its Passive option from the controller after this round.
