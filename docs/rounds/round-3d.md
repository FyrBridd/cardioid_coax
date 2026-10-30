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
the woofers (item 9), and since the review it decides on the compensator at the same time as the crossover points.
In this run the passive network meets every check listed in item 9. Every number below is from that run.

## What was run

- **Code:** branch `FyrBridd/passive-xo`, commit `dc9d061`. The BEM fields came from the cache, so there were no new
  solves. Outputs are not committed.
- **Command:** `uv run speaker passive --params params/round-3.toml --out out/passive/round-3-design --trusted-limit 2851`.
  It took 4 min 35 s (the search 256 s: the grid now fits every pair twice) and exited 0. The log is
  `out/passive/run.log`. `--trusted-limit` is now required: pass the trusted limit of the directivity run.
- **Outputs** (`out/passive/round-3-design/`): `report.html`, `network.txt`, `parts.csv`, `semi_active.txt`,
  `scores.json`, `schematic_bass.svg`, `schematic_coax.svg`, `links.svg`, `vituixcad_netlist.txt` and `plots/`
  (now with `electrical.png`).

## What changed from the plan while building it

1. **The level is tied to the woofers'.** This is a soft target, not a hard rule. The fit's level is the weighted
   mean of the 23° response (1/3-octave smoothed) over 100 Hz to 10 kHz: weight 1 from 100 Hz to the trusted limit
   (2851 Hz), and 0.25 from there to 10 kHz. The fit is penalised when that level falls more than 1 dB under the
   woofers' own level (80-150 Hz). **In this run it sits 1.6 dB under**, 0.6 dB past the 1 dB it may drop for free.
   The sensitivity is measured differently (23°, 300-1000 Hz), and is 1.5 dB under what the woofers give alone. A
   first preview run (not in this round's outputs) with the level left free threw away about 7 dB through a 51 Ω
   series resistor.
2. **The woofers must stay quiet under the coax LF.** From twice the bass crossover (340 Hz) up to 600 Hz they
   should sit at least 10 dB under the coax LF. This is a penalty in the fit, not a guarantee. Without it, the passive
   woofers ran about 8 dB under the coax out to about 800 Hz (an earlier preview run). That filled the coax LF's
   low-end droop (the cardioid chamber leaves the coax LF about 9 dB down at 150 Hz) and let non-cardioid sound into
   the band the cardioid section must own. **The delivered network nearly meets it, and misses at the bottom of the
   band.** Recomputed from this run (23°, 1/3-octave smoothed, as the fit does): the woofers sit 9.7 dB under the
   coax LF at 356 Hz, so 0.3 dB short, 1 of 11 grid points; 11.2 dB at 378 Hz. From 400 Hz up they are 12.8 dB under
   or more, reaching 21.4 dB at 600 Hz.
3. **A coax LF contour (shelf) network is offered.** It is R_lc with L_lc across it, before the coax LF pad. The run
   before the bass-peak guards used it (5.6 Ω with 0.7 mH, flat23 1.91 to 1.75 dB). **This run leaves it out:** it
   made flat23 worse (1.51 to 1.68 dB).
4. **Semi-active: the coax DSP channel gets EQ.** The design said gain only. The EQ is fitted before the DSP
   high-pass, as the active fit does, and only inside 150 Hz to 2 times the LF-HF point, cut at the trusted limit
   (150 Hz to 2400 Hz).
5. **Terminal links K1 to K5.** Each network is isolated at its boundary, so the same boards serve all three modes
   (`links.svg`, and the table at the end of `network.txt`). **K3, K4 and K5 are two-pole:** each breaks both the hot
   and the return between a network and its driver. In the active modes every terminal pair then reaches its driver
   alone. Without that, the drivers' returns stayed joined at the network ground, and in passive mode the inverted
   woofers have their + on that ground: wiring the BASS terminals "straight to the woofers" could have shorted an
   amplifier channel. The schematics now draw each driver last on its node, behind its link (drawn as a switch
   marked "K3 (2-pole)" and so on). So the Zobel and the compensator are visibly on the network side of K3, and the
   HF's RLC shunt on the network side of K5. **Wire them that way.** If the compensator ended up on the woofers'
   side of K3, the bass amplifier in the active modes would see 2.6 Ω at about 30 Hz.
6. **Impedance and tolerance rules.** The floor is 3.5 Ω. The tolerance check draws every inductor and capacitor at
   plus or minus 5 % (200 draws). If it fails, the fit is redone with the floor 0.1 Ω higher (up to 3.8 Ω). In this
   run it passed at 3.5 Ω. The series-woofer fallback is tried only if no parallel network clears 3.5 Ω.
7. **New: an electrical stopband rule.** The first network (L_b1 8 mH, C_b1 3.3 µF) had a bass transfer of
   **+3.0 dB at 2.43 kHz**: C_b1 was far too small to shunt the woofers, and it resonated with their voice-coil
   inductance. The fit could not see it. The bass field is exactly zero above the trusted limit, the bass-under term
   stops at 600 Hz, and the flatness is 1/3-octave smoothed. Now each branch's voltage transfer, on the modelled
   driver loads and independent of the acoustic data, must stay 12 dB under its level at its crossover point, from
   4 times past that point on: here the bass from 4 × 170 = 680 Hz up, the coax LF from 4.9 kHz up and below 42 Hz,
   and the coax HF below 300 Hz. The fit penalises every dB over that at 20 per dB (`STOP_WEIGHT`). On this design
   (before item 8) weights 3 and 10 left the bass 0.7 and 1.8 dB over the line near 4 times its crossover point.
   Weights 15 to 100 held it within 0.3 dB. 20 is one step above the lowest that held. `plots/electrical.png` shows
   each branch with the limits.

   **What the rule cost when it was added.** Flatness at 23° went from 1.43 dB (first network) to 1.75 dB. The
   bass-to-coax point moved from 150 Hz to 170 Hz. The level dropped by about 2 dB. It did not force the series
   fallback. **In this run every stopband holds, two of them narrowly:**
   - the bass is −13.9 dB re its 170 Hz level at 686 Hz (limit −12 dB);
   - the coax HF is −13.1 dB re its 1200 Hz level at 298 Hz;
   - the coax LF is −16.1 dB below 42 Hz and −25.0 dB from 4.9 kHz.
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

   Both are penalties, 1000 per dB at each frequency point (`PEAK_WEIGHT`). Before the compensator existed I tried
   weights from 20 to 10000 on the whole search. No weight met the +1 dB cap, and the result moved a lot with the
   weight (500 gave flat23 6.06 dB, 2000 gave 1.99 dB). The best, at 1000, held the cap only by making the bass
   low-pass a gentle slope into a shelf. That missed the stopband rule of item 7 by 6.3 dB, and whole-system rear
   rejection fell to 7.93 dB. The owner chose to fix that with the compensator (item 9).

   **Guards not covered: a resonance behind a pad.** The gain cap reads each driver's terminal voltage. A resonance
   at a node inside a network, ahead of an L-pad, does not show there. In this run the coax high-pass (C_h1, L_h1)
   peaks at node cx, ahead of the coax LF's low-pass and pad: **+6.0 dB at 116 Hz** (+5.0 dB at 108 Hz in the first
   compensator run, commit d1a4a7b). The pad brings the coax LF down to −7.1 dB at 108 Hz and −7.0 dB at most
   (112 Hz), so the cap does not see it. The coax LF still gets more drive near 110 Hz than anywhere else in its band
   (see Open, excursion).
9. **New, the owner's choice: an impedance-peak compensator across the woofers.** The woofer pair in its sealed box
   has an impedance peak of about **29 Ω at 40 Hz** (the box resonance), and falls to 3.26 Ω at 147 Hz. Above the
   peak the woofers look like a capacitor, and that is what the bass coil resonated with (item 8). The compensator is
   the standard fix: a resistor, a coil and a capacitor in series (R_bc, L_bc, C_bc), wired across the woofers on the
   network side of K3. Near the box resonance it takes current the woofers do not, so the network sees a much flatter
   load and a normal low-pass can work into it. Its textbook values come from the woofers' Re, fc, Qec and Qmc:
   3.6 Ω, 9.2 mH and 1600 µF here. With those, the pair's 20-200 Hz impedance would run from 2.7 to 3.3 Ω instead of
   3.3 to 29 Ω. That is under the 3.5 Ω floor, so the fit detunes it: **2.7 Ω, 34.5 mH and 820 µF**, tuned to 30 Hz.
   With these the woofer load runs from 2.7 Ω at 30 Hz to 8.6 Ω at 53 Hz (31 Ω without). The series coil lifts the
   amplifier's load to 3.57 Ω at 27 Hz.

   **How the search decides it.** The grid fits every pair of crossover points twice, once without and once with
   the compensator, and picks among all 60 candidates with the usual rules. A first version offered the compensator
   only after the points had been picked without it. That run kept it, but at points chosen for the bare network, and
   a review rerun showed the result swung with the guard weight. With the grid deciding, the pick was **170 / 1200 Hz
   with the compensator** (flat23 1.51 dB at the grid). The options then kept the Zobel (1.51 to 1.35 dB) and the RLC
   shunt across the coax HF (1.35 to 0.64 dB). They left out the contour, a 3rd-order bass, a 3rd-order coax LF and
   the notch.

   **Stability check.** I ran the whole search at guard weights 500, 1000 and 2000:
   - **1000 and 2000 gave the same network** (170 / 1200 Hz, compensator, Zobel, HF shunt), and it met every check.
   - **500 chose 150 / 1400 Hz** with the compensator but no other option. It met every check except rear
     rejection: whole_rear 7.90 dB (limit 9.5). Its flat23 was 1.35 dB and presence lift 4.29 dB.

   So the result is stable from 1000 up and not below. I kept 1000. Before this change, with the compensator
   offered only after the grid, the reviewer's runs at weights 100, 500 and 2000 all dropped it, and none of them
   met every check.

   **The result** (passive, weight 1000):
   - Impedance: lowest 3.51 Ω (at 20 kHz), and the tolerance check passed (worst draw 3.34 Ω).
   - Bass stopband: −13.9 dB at 686 Hz.
   - The woofers get at most **+1.0 dB, at 64 Hz** (cap +1 dB).
   - The sum below 100 Hz rises at most **0.5 dB above the fit's level**, at 67 Hz (limit 1 dB): no hump.
   - whole_rear **9.72 dB**, and flat23 **0.64 dB**, the flattest of the three modes.

   **What it costs.** Three more parts, two of them large:
   - L_bc is a 34.5 mH cored coil. It must carry 7.8 A peak without saturating, and its DCR heats by 25 W at full
     voltage.
   - C_bc is an 820 µF bipolar electrolytic: 56 V peak (a 100 V part) and **5.5 A rms ripple** at full voltage.
   - R_bc (2.7 Ω) dissipates 9.1 W with a third of the amplifier's voltage and 82 W at full voltage, both at
     about 28 Hz.

   The Zobel (1.8 Ω with 150 µF) is a second shunt across the woofers. Sensitivity is 79.3 dB, 1.5 dB under the
   woofers alone.

## The passive network

**Crossover points: bass to coax 170 Hz, coax LF to HF 1200 Hz.** These were picked from 60 grid candidates: 30
pairs, each with and without the compensator. The network uses the woofers in parallel, the bass inverted and the
coax HF inverted. One amplifier, one pair of terminals (IN), and links K1 to K5 all closed. Every value is a
standard value (capacitors E12, resistors E24).

| Part | Value | Type and rating | DCR | At full voltage | 1/3-voltage power |
|---|---|---|---|---|---|
| L_b1 | 7.5 mH | cored inductor | 0.38 Ω | 10.1 A peak (saturation at least this); DCR heat 19.6 W | |
| C_b1 | 22 µF | film capacitor, 63 V | | 39.5 V peak, 0.26 A rms | |
| R_bz | 1.8 Ω | resistor (Zobel), 5 W | | 5.8 W | 0.64 W |
| C_bz | 150 µF | bipolar electrolytic (Zobel), 63 V | | 39.3 V peak, 1.79 A rms | |
| R_bc | 2.7 Ω | resistor (compensator), listed 10 W: **buy 20-25 W** | | 82 W | 9.09 W |
| L_bc | 34.5 mH | cored inductor (compensator) | 0.82 Ω | 7.8 A peak (saturation at least this); DCR heat 24.9 W | |
| C_bc | 820 µF | bipolar electrolytic (compensator), 100 V | | 56.4 V peak, **5.50 A rms** | |
| C_h1 | 100 µF | film capacitor, 160 V | | 75.4 V peak, 7.14 A rms | |
| L_h1 | 20 mH | air-core inductor | 2.24 Ω | 4.9 A peak; DCR heat 27.1 W | |
| L_l1 | 4 mH | air-core inductor | 1.00 Ω | 2.2 A peak; DCR heat 2.4 W | |
| C_l1 | 5.6 µF | film capacitor, 160 V | | 69.9 V peak, 1.35 A rms | |
| R_ls | 33 Ω | resistor, 10 W | | 46 W | 5.10 W |
| R_lp | 24 Ω | resistor, 5 W | | 5.2 W | 0.58 W |
| C_t1 | 3.9 µF | film capacitor, 160 V | | 70.8 V peak, 7.25 A rms | |
| L_t1 | 1.5 mH | air-core inductor | 0.61 Ω | 0.34 A peak | |
| C_t2 | 180 µF | film capacitor, 50 V | | 0.8 V peak, 7.24 A rms | |
| R_ts | 1 Ω | resistor, 10 W | | 52 W | 5.82 W (semi-active) |
| R_tp | 2.4 Ω | resistor, listed 10 W: **buy 20-25 W** | | 89 W | 9.92 W (semi-active) |
| R_tr | 3 Ω | resistor (HF shunt), 5 W | | 2.6 W | 0.29 W |
| L_tr | 0.2 mH | air-core inductor (HF shunt) | 0.22 Ω | 1.3 A peak | |
| C_tr | 12 µF | film capacitor (HF shunt), 50 V | | 4.3 V peak, 0.93 A rms | |

Full-voltage figures are a sine with the amplifier's full 25 V rms on the input, the largest from 20 Hz to 20 kHz,
in either mode. A capacitor's rating is the first of 50, 63, 100, 160 and 250 V at or above 1.5 times its peak
voltage. The ripple current matters for the electrolytics: check C_bc's 5.5 A rms and C_bz's 1.8 A rms against the
datasheets. Parts per pair are in `parts.csv` (columns `peak_a`, `peak_v`, `ripple_a`).

**What each part does.**
- **Bass (L_b1, C_b1; the Zobel R_bz, C_bz; the compensator R_bc, L_bc, C_bc).** A 2nd-order low-pass (coil, then
  capacitor to ground) into the woofers, with two shunts across them on the network side of K3. Its transfer
  (`plots/electrical.png`):
  - −2.5 dB at 30 Hz and at most **+1.0 dB at 64 Hz**;
  - −3.5 dB at 100 Hz and −9.0 dB at 170 Hz;
  - then falling: −13 dB at 300 Hz, −23 dB at 680 Hz, −28 dB at 1 kHz, −36 dB at 2 kHz and −40 dB at 3 kHz.

  It meets the stopband rule (−13.9 dB re 170 Hz at 686 Hz, limit −12 dB).
- **Coax high-pass (C_h1 from IN to node cx, L_h1 from cx to ground).** It feeds both coax networks, the LF and the
  HF. It resonates: node cx is **+6.0 dB over the input at 116 Hz** and +4.4 dB at 100 Hz (item 8, "a resonance
  behind a pad").
- **Coax LF.** A low-pass (L_l1, C_l1) and a pad (R_ls, R_lp), no contour. The coax LF gets its most voltage,
  **−7.0 dB, at 112 Hz** (the cx resonance through the pad). It gets −7.9 dB at 100 Hz, −10.2 dB at 150 Hz, and
  −12 dB around 1 kHz.
- **Coax HF.** A high-pass (C_t1, L_t1, C_t2), a pad (R_ts, R_tp) and a series RLC shunt across the driver (R_tr,
  L_tr, C_tr, tuned near 3.2 kHz). The shunt takes the edge off the HF's 2-3.3 kHz hump: flat23 went from 1.35 to
  0.64 dB with it. C_t2 (180 µF) is large enough to be nearly a wire, so the high-pass is close to 2nd order. A
  180 µF film capacitor is a big, costly part.

**What the search kept and left out:**
- The grid picked the compensator with the points (above).
- Kept the Zobel across the woofers (flat23 1.51 to 1.35 dB) and the RLC shunt across the coax HF (1.35 to
  0.64 dB).
- Left out the coax LF contour (1.51 to 1.68 dB), a 3rd-order bass low-pass (1.35 to 1.52 dB), a 3rd-order coax LF
  low-pass (1.35 to 1.37 dB) and the coax LF notch (1.35 to 1.39 dB): each was no better.

**Impedance and level:**
- Minimum impedance is **3.51 Ω at 20 kHz** (+0.01 Ω over the floor). In the bass band (20-500 Hz) it is **3.57 Ω
  at 27 Hz**. `plots/impedance.png` shows the 3.5 Ω line. The floor is the SA550's minimum.
- **Semi-active impedance.** The coax channel's minimum is **3.46 Ω at 20 kHz, just under 3.5 Ω.** The bass channel
  (bare woofers in parallel; the compensator is behind K3, so it is out) dips to **3.26 Ω at 149 Hz**. Neither
  channel is covered by the fit or the tolerance check, because the amplifiers for semi-active are not chosen. Pick
  ones that are happy at about 3 Ω.
- Sensitivity is **79.3 dB at 2 m for 2.83 V** (23°, 300-1000 Hz). That is 1.5 dB below what the woofers give alone.
- **Tolerance check passed at 3.5 Ω.** The worst of 200 draws was 3.34 Ω (limit 3.3 Ω). Flatness (flat23) was
  0.64 dB nominal and 0.91 dB at the 95th percentile (limit +0.5 dB). Whole_rear ran 9.50 to 9.95 dB (5th to 95th
  percentile).
- No series-woofer fallback was needed: a parallel network cleared the floor.

**Power.** Resistor ratings follow the rule in `network.txt`: the first of 5, 10, 15, 20, 30, 50 W at or above the
dissipation with one third of the amplifier's voltage (roughly what music averages). Where that figure reaches 85 %
of the rating, `network.txt` and `parts.csv` now say to buy 2 to 2.5 times the rating:
- **R_bc** (2.7 Ω): 9.09 W against 10 W, so **buy 20-25 W**. Its worst case is at about 28 Hz, where a sustained bass
  note is close to a sine, so the sine figure is not conservative there. At full voltage it would take 82 W.
- **R_tp** (2.4 Ω): 9.92 W against 10 W in semi-active mode (9.69 W passive), so **buy 20-25 W**. Its worst case is
  at 20 kHz, where music has far less energy, so the sine figure is conservative for R_tp.
- R_ts at 5.82 W and R_ls at 5.10 W, both against 10 W. Every other resistor is under 1 W.

At full voltage R_tp would take 89 W, R_bc 82 W, R_ts 52 W and R_ls 46 W.

**Coils and capacitors.** Two cored coils: L_b1 (7.5 mH, 10.1 A peak) and L_bc (34.5 mH, 7.8 A). Each must be
rated for at least that saturation current, and each heats by 20 to 25 W in its DCR at full voltage. L_h1 (20 mH
air core, 2.24 Ω, 27 W of DCR heat at full voltage) is large. C_bc (820 µF) and C_bz (150 µF) are bipolar
electrolytics. C_h1 (100 µF) and C_t2 (180 µF) are large film capacitors.

### Warnings in network.txt

- **L_h1 is at the 20 mH limit.** The fit wanted an even larger coil for the coax high-pass, but 20 mH air-core is
  already 2.24 Ω of DCR and would be large and costly. The response near this part is a compromise.
- No line is marked **OVER**. The largest branch gain is +1.0 dB at 65 Hz, on its +1 dB limit.

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
| bass | −2.62 dB | 1.65 ms | inverted | 80 Hz +1.18 dB Q1.65; 214 Hz −1.77 dB Q1.22 |
| coax | 0 dB | 0 | normal | 150 Hz **+12.0 dB** Q0.75 (at the fitter's limit); 205 Hz −3.84 dB Q0.92; 588 Hz −2.28 dB Q1.55; 1734 Hz −1.36 dB Q1.17 |

- **The 150 Hz PEQ is at the fitter's limit again.** These coax boards pass the coax LF only −10.8 dB at 150 Hz in
  semi-active mode, so the DSP has a lot of the cardioid chamber's low-end droop to fill.
- **Stacked boost.** The four coax PEQs together (my sum of their curves) peak at **9.2 dB at 140 Hz**. That needs
  the usual care: amplifier clipping and the coax LF's excursion near 150 Hz should be checked (see Open). Normalise
  in the DSP so the largest boost does not clip, as the file says.

## Three modes, trusted scores

Trusted band (150 Hz to 2851 Hz; bass DI from 120 Hz). Lower is better except the rear and cardioid rows. The first
three columns are earlier passive networks, for comparison: the first network, the one with the stopband rule but
before the bass-peak guards, and the one with the guards but no compensator. Semi-active and active are from this
run.

| Score (dB) | Passive, first network | Passive, before the guards | Passive, guards alone | Passive | Semi-active | Active |
|---|---|---|---|---|---|---|
| 23° flatness | 1.43 | 1.75 | 1.42 | **0.64** | 1.24 | 1.16 |
| Window spread | 3.60 | 3.73 | 3.99 | 3.45 | 3.45 | **3.40** |
| Presence lift, 2-2.85 kHz | 4.22 | 4.04 | 4.32 | **3.66** | 4.04 | 3.90 |
| DI smoothness, bass, 120-500 Hz | 0.94 | 0.69 | 1.32 | **0.68** | 1.42 | 1.66 |
| DI smoothness, LF-HF, 500-2851 Hz | 0.77 | 1.11 | 1.20 | **0.95** | 1.09 | 1.08 |
| Whole-system rear rejection, 150-300 Hz (higher is better) | 10.44 | 11.62 | 7.93 | 9.72 | **9.98** | 9.88 |
| Cardioid hard / aim (coax section) | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 | 17.15 / 12.36 |

Bold marks the best of this run's three modes. The coax section's cardioid scores are the same because the coax is
the same box. Semi-active's scores moved with the coax boards. With the guards alone it had flat23 0.81, window
spread 3.97, presence lift 4.23, DI LF-HF 1.12 and whole_rear 9.87 dB.

## Where passive falls short, and why

- **Flatness: passive now leads.** Passive is 0.52 dB flatter than active at 23° (0.64 against 1.16 dB) and 0.60 dB
  flatter than semi-active (1.24 dB). The coax HF shunt does most of that. The owner can judge all three from
  `plots/crossover_23_*.png`.
- **The 2 to 5 kHz window: passive is now the lowest of all.** Its presence lift is **3.66 dB**:
  - 0.24 dB under active (3.90 dB) and 0.38 dB under semi-active (4.04 dB);
  - under every earlier passive network (4.22, 4.04 and 4.32 dB).

  It is still well above the 0.5 dB the score calls a pass. The lift is the coax's, and was already there in
  round 3c. Window spread is 3.45 dB in passive and semi-active, against 3.40 dB in active. Above 2851 Hz the scores
  are low confidence and are not ranked. The +4 dB at 3.3-3.8 kHz in `plots/crossover_23_passive.png` lies in the
  low-confidence cabinet-to-puck blend, and active shows it too.
- **Below the fitted band: no hump.** Below 100 Hz the passive sum (23°, 1/6-octave smoothed) is 76.8 dB at 50 Hz,
  78.9 dB at 60 Hz, 79.7 dB at 67 Hz and 77.3 dB at 100 Hz, then 79.0 to 79.4 dB from 150 to 300 Hz. Against the
  fit's level (1/3-octave, as the fit sees it) the highest point below 100 Hz is 0.5 dB above, at 67 Hz.
- **Rear rejection.** The whole_rear score is the minimum over 150 to 300 Hz. Passive is 9.72 dB against 9.98 and
  9.88 dB. Across the band it is lower: at about 240 Hz the passive rejects about 13.2 dB, and semi-active and active
  about 18.5 dB (`plots/whole_rear_compare.png`). Below 150 Hz the passive rejects more.
- **Sensitivity.** 79.3 dB is 1.5 dB below the woofers alone. The active figure is not comparable, as each way has
  its own amplifier gain.
- **Size and cost.** The passive boards need:
  - two cored coils (7.5 mH at 10 A, 34.5 mH at 7.8 A) and a 20 mH air-core coil;
  - an 820 µF and a 150 µF electrolytic, and a 100 µF and a 180 µF film capacitor;
  - two resistors to buy at 20-25 W.

  None of these is small or cheap. Parts cost was not estimated.

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

- **Whether the owner wants passive at all.** In this run passive is the flattest of the three modes, and has the
  least 2-5 kHz lift. It needs one amplifier and no DSP. Against that, it has a little less rear rejection through
  the middle of the cardioid band, and the parts above. Semi-active needs two channels and a DSP on the coax only.
- **The rule conflict is resolved by the compensator, from guard weight 1000 up** (item 9). At 1000 and 2000 the
  network meets every check. At 500 it picks other points and misses rear rejection. The compensator's values depend
  on the woofers' modelled resonance, so measured driver data may move the result.
- **Excursion in passive mode.**
  - The woofers get at most +1.0 dB, at 64 Hz.
  - The coax LF gets its most voltage near 110 Hz, from the cx resonance: −7.0 dB at 112 Hz, −7.9 dB at 100 Hz.
    No guard limits that (item 8).
  - Judging excursion needs the drivers' Xmax and an excursion model. The coax LF near 110 Hz is the one to check.
- The stacked 9.2 dB semi-active coax boost at 140 Hz, with its 150 Hz PEQ at the fitter's limit, needs a check of
  the coax LF's excursion and amplifier headroom.
- **Buy these parts with margin.**
  - R_bc and R_tp at 20-25 W.
  - L_b1 and L_bc rated for 10.1 A and 7.8 A.
  - C_bc for 5.5 A rms of ripple at 100 V.
- Everything above depends on unmeasured drivers. Rerun after the coax and woofers are measured.
- The VituixCAD netlist (`vituixcad_netlist.txt`) has not been imported or checked by the owner.
- The owner's report page gets its Passive option from the controller after this round.
