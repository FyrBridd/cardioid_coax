# Round 1: cardioid first, then the front shape

Date: 2026-10-27/28. The first optimisation round (spec §2). Stage 1 searched
the vent position and slot count, tuning the felt and the coax chamber volume
for each. Stage 2 tried front shapes on the stage 1 winner. The best two went
to full-resolution finals.

## What was run

- **Code commit:** `93a3d55` (branch `FyrBridd/plan2-search-review`).
- **Round file:** `rounds/round-1.toml`, sha256 `819eb9cf…`.
- **BEM code version** (cache key): `d00d81a6…`.
- **Command:** `uv run speaker round rounds/round-1.toml`, then the mesh
  check and full rerun of the recommendation (below).
- **Outputs** (not committed): `out/rounds/1/report.html`, with
  `notes.md`, `results.json`, `finals.json`, each variant's `params.toml`,
  `params_tuned.toml` and `scores.json`, and the finals' own reports.

| Stage | Settings | Ranked by | Variants | Rejected | Errors | Time |
|---|---|---|---|---|---|---|
| vents | vents (to 600 Hz) | tuned cardioid margin | 15 + stage base + round 0 | 0 | 0 | 1 h 45 min (est. 1 h 45 min) |
| shape | screen (5 per λ, to 2715 Hz) | directivity, behind the cardioid gate | 14 + stage base + round 0 | 0 | 0 | 5 h 09 min (est. 6 h 09 min) |
| finals | full run (6 per λ, with the puck) | – | 2 | – | 0 | 4 h 01 min (est. 3 h 50 min) |
| Total | | | | | | 10 h 55 min (14:24–01:19) |

The six shape rows the plan listed after `w330-fa25-ramp-m240` (thicker outer
baffle, cone flare, coax offset in and out, 60 mm top margin, lower bass
heights) were commented out of the round file to fit the 12 h budget. They
were not run.

## The cardioid now passes

Moving the vents forward fixes round 0's main failure. The tune sets the felt
resistance and the coax chamber net volume for each vent layout, for the
worst case over the ±30 % felt spread.

| Design | Vents (from front, slots) | R (rayl) | V (L) | Hard, 150–300 Hz | Aim, 300–500 Hz |
|---|---|---|---|---|---|
| Round 0 | 270 mm, 12 | 190 | 10.0 | 5.3 dB **fail** | 5.1 dB warn |
| Round 0 front, tuned (stage base) | 270 mm, 12 | 236 | 14.5 | 7.7 dB **fail** | 6.7 dB |
| Vents winner | 120 mm, 9 | 115 | 16.0 | **16.6 dB** pass | 11.5 dB |
| Recommendation (finals) | 120 mm, 9 | 115 | 16.5 | **16.3 dB** pass | 11.6 dB |

- 12 of the 15 vent layouts pass the 10 dB gate. The three that fail are the
  furthest back (240 mm with 12 slots, 280 mm with 9 or 12).
- Position matters most: taking the best slot count at each position, each
  40 mm forward gains 1.2–1.9 dB (10.0, 11.8, 13.6, 14.8, 16.6 dB from
  280 mm to 120 mm). The slot count moves the result by 0.6–2.8 dB at a
  given position; 6 slots is best from 160 mm back, 9 at 120 mm.
- The rejection is still rising at 120 mm, the front edge of the grid.
- **Felt spread** (recommendation): hard / aim = 13.5 / 12.2 dB at R × 0.7,
  16.3 / 11.6 dB at R × 1, 12.9 / 11.0 dB at R × 1.3. The gate holds across
  the spread.
- The worst point is still 300 Hz, straight behind. At 50 m the hard score
  is 15.1 dB, so it is not a near-field effect.
- The rear rejection curve has a deep null near 190 Hz (about 42 dB). The
  felt spread moves it, but the band minimum stays at 12.9 dB or more.
- The lumped model still disagrees with the 3D result by up to 37.6 dB in
  150–500 Hz. The ranking uses the 3D result only.
- **Side effects of the larger chamber:** the printed fillers shrink from
  8.8 L to 3.2 L (16.5 L net from a 25.2 L gross chamber). The bass section
  is unchanged at 54.1 L net (56.4 L at 330 mm wide).

## The front shape

All 14 shape variants (plus the stage base) carry the vents winner, so they
all pass the gate (16.3–16.7 dB). They are ranked on directivity: window
spread plus presence lift in 0.25 dB buckets, then DI smoothness, then 23°
flatness. Screen scores, trusted, in dB:

| Rank | Variant | Window spread | Presence lift | Sum | DI (bass + LF–HF) | flat23 |
|---|---|---|---|---|---|---|
| 1 | w330-sq-m290 | 3.29 | 3.50 | 6.79 | 2.07 | 3.05 |
| 2 | stage base (320, square, m290) | 3.56 | 3.73 | 7.30 | 2.08 | 2.98 |
| 3 | w330-ro20-m280 | 3.64 | 3.72 | 7.36 | 2.20 | 3.07 |
| 4 | w330-fa15-step-m260 | 3.53 | 3.93 | 7.46 | 2.45 | 3.17 |
| 5 | w320-ro15-m280 | 3.75 | 3.87 | 7.62 | 2.21 | 3.01 |
| 6 | w330-fa15-ramp-m240 | 3.46 | 4.20 | 7.66 | 2.86 | 3.33 |
| 7 | w330-ch20-m280 | 3.90 | 4.02 | 7.93 | 2.24 | 3.06 |
| 8 | w330-ro30-m260 | 3.85 | 4.07 | 7.92 | 2.44 | 3.15 |
| 9 | w330-fa25-step-m260 | 3.79 | 4.18 | 7.97 | 2.48 | 3.16 |
| 10 | w320-ro35-m240 | 4.07 | 4.25 | 8.32 | 2.83 | 3.27 |
| 11 | w330-ro40-m240 | 4.04 | 4.39 | 8.44 | 2.91 | 3.31 |
| 12 | w330-fa25-ramp-m240 | 3.85 | 4.53 | 8.38 | 2.93 | 3.32 |
| 13 | w330-ch30-m260 | 4.32 | 4.56 | 8.88 | 2.51 | 3.12 |
| 14 | w320-ch35-m240 | 4.57 | 4.78 | 9.35 | 3.00 | 3.19 |
| 15 | w330-ch40-m240 | 4.76 | 5.00 | 9.75 | 2.99 | 3.27 |
| – | round 0 (reference) | 3.56 | 3.75 | 7.31 | 2.11 | 3.02 |

What this shows:

- **No edge treatment beat the square edge.** The bigger the roundover,
  chamfer or facet, the worse the window. Chamfers did worst.
- **But edge and mouth are confounded.** Every edge variant also had a
  smaller flare mouth (240–280 mm against 290 mm); the round file paired
  them. This round cannot say which of the two made the window worse.
- **The vent move did not change directivity.** Round 0 and the stage base
  differ by at most 0.04 dB on every directivity score.

## Finals and the recommendation

The finals rerun the top two at 6 elements per wavelength, with the puck and
the full 16 kHz grid. Trusted scores, in dB:

| Score | w330-sq-m290 (rank 1) | stage base (rank 2) | Round 0b (earlier code) |
|---|---|---|---|
| cardioid_hard | 16.27 pass | 16.63 pass | 5.28 fail |
| cardioid_aim | 11.56 pass | 11.47 pass | 5.13 warn |
| window_spread (150–2851) | 3.55 warn | 3.80 warn | 4.05 warn |
| presence_lift (2000–2851) | 4.19 warn | 3.91 warn | 5.26 warn |
| flat23 (150–2851) | 3.60 | 3.71 | 4.14 |
| di_smooth_bass | 0.62 | 0.63 | 0.54 |
| di_smooth_lfhf | 1.44 | 1.44 | 1.60 |
| whole_rear | 3.58 | 3.45 | 2.71 |

**Recommendation: w330-sq-m290**, the rank 1 final. It is round 0 with the
vents at 120 mm and 9 slots, R ≈ 115 rayl, V = 16.5 L, and the width at
330 mm. Square edges, 290 mm elliptical flare, everything else unchanged.

The width is a free choice, not a finding:

- On the screen, 330 mm was 0.5 dB better than 320 mm on window spread plus
  presence lift. On the finals 320 mm is 0.04 dB better on that sum. Both
  are in the same 0.25 dB bucket, and 330 mm wins the DI tie-break by
  0.004 dB.
- Every score is within 0.4 dB between the two. 330 mm has 0.25 dB less
  window spread; 320 mm has 0.28 dB less presence lift and 0.36 dB more
  cardioid margin.
- If you would rather keep the 320 mm drawings, the stage base is the same
  design and equally good.

**Directivity against round 0b is not like-for-like.** Window spread and
presence lift look 0.5–1.1 dB better than round 0b, but the rescored round 0
on the screen matches the stage base. The gain comes from the Task 2 changes
to the coax de-embedding (the recess ratio R and the HF reference), not from
the design. The real picture is unchanged from round 0: the 0–10° directions
rise about 3–3.5 dB over 23° by 2–2.8 kHz, and the window stays within
±1 dB only up to about 710 Hz (round 0b: 898 Hz, earlier code).

**Crossover gains** (fixed LR4 at 170 / 1400 Hz, matched at 23°): bass
+11.6 dB, coax LF 0, coax HF −4.9 dB. Round 0 needed +14.2 dB on the bass
and −2.3 dB on the HF.

## Trusted limit of the recommendation

**By the rule: 177 Hz.** The mesh check of w330-sq-m290 (6 against 8
elements per wavelength, as in round 0) failed the 198 Hz band by 2.41 dB,
so the limit fell to that band's lower edge. The full rerun then refused to
start (it needs a limit of at least 500 Hz). The round's finals report,
with the inherited 2851 Hz, stays as the recommendation's report.

**Why: one point in a very deep null.** I recomputed the compared curves
from the cached fields (a scratch script, not committed):

- The whole 2.41 dB is the coax-section rear rejection at ±165° and
  178 Hz: 58.1 dB on the run mesh, 55.7 dB on the fine mesh. The 0.28 dB in
  the 157 Hz band is the edge of the same null.
- At 56–58 dB the rear pressure is about 1/800 of the front, so the two
  meshes differ by a negligible absolute amount. The rule compares the
  rejection in dB, which amplifies it. Round 0 had no null this deep, so
  the rule never met one.
- Every other compared curve (23°, listening window, early reflections,
  sound power, H 0–90°) agrees within 0.07 dB up to 2851 Hz (worst at
  2136 Hz).
- With the rejection capped at 30 dB before the comparison, the worst
  difference is 0.09 dB and the limit is **2851 Hz**, the same as round 0.
  A 40 dB cap gives the same limit (worst 0.29 dB).
- The check again stopped at the 4 kHz band on memory (about 24.8 GB
  against 22 GB).

So the scores stand. The cardioid scores use the band minimum (16.3 dB at
300 Hz), far from the null. Formally, though, this design's own limit is
177 Hz until the rule caps the rejection before comparing (see Still open).

## The recess finding (round 0 open question 1)

Task 2 measured the coax cone's own recess: the BEM on-axis pressure of the
cone in its frame recess over that of a flat piston, on a 1 m disc baffle.
That ratio R now divides the datasheet before the 3D model adds the recess
back, so it is no longer counted twice.

| Frequency | R |
|---|---|
| 300 Hz | +1.1 dB |
| 570 Hz | +3.7 dB |
| 756 Hz (peak) | +4.3 dB |
| 1 kHz | +2.8 dB |
| 1.9 kHz | −11.3 dB |
| 3 kHz | −0.2 dB |

- **Open question 1 (the 550–650 Hz peak).** Round 0's cone field was
  +13.5 dB over a point source at 570 Hz, about 7.5 dB more than a flat
  baffle explains. R accounts for about 3.7 dB of that, so about half the
  excess was the double count and about half remains. Read off the new
  spinorama, the on-axis peak is about 94.5 dB at 2 m, against about 98 dB
  in round 0. The rest is either real (the recess and flare acting as a
  short horn) or the guessed cone profile. Measuring the cone settles it.
- **Open question 2 (the 1.8–1.9 kHz notch).** R has a −11.3 dB notch at
  the same frequency, so the notch is the recess itself. With R divided out,
  the result there now follows the datasheet. It is above the 1.4 kHz
  crossover in any case.
- The smaller coax LF peak probably explains most of the 2.6 dB drop in the
  bass gain.

## Still open

1. **Vents further forward.** 120 mm was the edge of the grid and the
   rejection was still rising. Round 2 should try 80–100 mm, if the chamber
   layout allows.
2. **Edge against mouth.** Separate the two: a 20 mm roundover with the
   290 mm mouth, and the square edge with a smaller mouth.
3. **The six dropped rows** (outer baffle 30 mm, cone flare, coax offset
   ±10 mm, top margin 60 mm, lower bass heights), on the square-edged front
   rather than the 30 mm roundover the round file names.
4. **The 2–2.8 kHz window rise.** The main directivity problem left. No
   baffle edge tried here helped. It may need the flare or the crossover
   (plan 3).
5. **Screen against finals.** The screen's 0.5 dB width advantage vanished
   in the finals. Treat screen differences under about 0.5 dB as ties.
6. **The mesh check's rejection rule.** Cap the rear rejection (at, say,
   30 dB) before comparing meshes, so a deep null cannot set the trusted
   limit. Your call; it is a code change, and the rerun is minutes because
   every field is cached.
7. **Round 0's open questions 3–5** (the rear null at ±120°, the low bass
   below the axis, the puck detail above the trusted limit) are unchanged.
   The provisional assumptions of round 0 (cone depth 65 mm, HF exit
   25.4 mm, cutout 230 mm, fixed crossover, basket-cylinder coax) still
   stand.
