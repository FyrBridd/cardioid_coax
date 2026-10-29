# Round 3b: round 3 re-scored with the fitted crossover

Date: 2026-10-29. Plan 3 replaced the fixed crossover (ideal LR4 at 170 Hz and 1.4 kHz, one gain per
way) with an active crossover fitted to each design. This round reruns round 3 through it. It also
makes a full run, a CAD build and a VituixCAD export of the owner's choice, bass260-580
(`params/round-3.toml`). Everything here is simulated; the drivers are not bought.

## What was run

- **Versions:** code commit `329132c` (branch `FyrBridd/plan3-crossover-exports`), BEM code version
  `d00d81a6…`, scoring version 2. The round file `rounds/round-3.toml` (sha256 `74fa4dd0…`) is
  unchanged. Outputs are not committed.
- **Full run:** `uv run speaker directivity --params params/round-3.toml --out out/plan3/round-3-design --trusted-limit 2851`
  took 21 s, all cached.
- **Round:** `uv run speaker round rounds/round-3.toml --out out/rounds/3b`. The shape stage took
  3 min, cached (10 variants, 0 errors). The finals took 44 min.
- **CAD:** `uv run speaker build --params params/round-3.toml --out out/plan3/round-3-build` took 4 s.
- **The one new BEM solve** was the bass280-600 final, which the new ranking brought in. It took
  2652 s, with memory in use peaking at about 17 of 30 GB.

## The fitted crossover of bass260-580

Points: **bass to coax 170 Hz, coax LF to HF 1000 Hz** (search ranges 150–250 Hz and 1.0–2.0 kHz).
The acoustic target is LR4 at each point, flat at 23°.

| Way | Gain | Delay | Polarity | PEQ (f0, gain, Q) |
|---|---|---|---|---|
| bass | +11.5 dB | 1.25 ms | inverted | 80 Hz +1.2 dB Q1.8; 213 Hz −1.8 dB Q1.3 |
| coax LF | 0 dB | 0 | normal | 183 Hz +6.0 dB Q2.5; 585 Hz −5.3 dB Q0.9; 1394 Hz −2.7 dB Q3.4; 1893 Hz +12.0 dB Q0.5 |
| coax HF | −4.8 dB | 0.01 ms | normal | 1885 Hz −4.7 dB Q1.0; 3997 Hz −6.6 dB Q2.9; 5620 Hz +1.9 dB Q5.4; 11831 Hz +2.0 dB Q3.2; 12833 Hz +12.0 dB Q0.65 |

The delays are well inside their grids (±3 ms and ±1 ms). Over 150–2851 Hz the sum at 23° stays
within −1.5 / +1.6 dB of its 200–1000 Hz level (91.9 dB at 2 m for 2.83 V), with no notches. The
fixed crossover gave −3.9 / +4.1 dB. `crossover.txt` is the DSP setup sheet, with biquads at 48
and 96 kHz.

**Passive-network warnings:**
1. EQ boost: the coax LF needs 10.9 dB at 2 kHz and the coax HF 13.7 dB at 12.1 kHz (passive: 6 dB).
2. The 1000 Hz point is within an octave of the 750 Hz HF resonance (the limit is 1500 Hz).
3. The coax HF filters at 4.0, 5.6, 11.8 and 12.8 kHz are above the trusted limit (2851 Hz).

## Scores against the fixed-crossover final

bass260-580, full run, in dB, trusted unless marked:

| Score | Fixed (round 3 final) | Fitted | Change |
|---|---|---|---|
| flat23 (150–2851 Hz) | 4.01 | **1.49** | −2.52 |
| window_spread (150–2851 Hz) | 3.48 warn | 3.36 warn | −0.13 |
| presence_lift (2000–2851 Hz) | 4.10 warn | 4.12 warn | +0.02 |
| presence_lift (2851–5000 Hz, low confidence) | 6.08 | 5.54 | −0.54 |
| di_smooth_bass | 0.55 | **1.35** | +0.80 (worse) |
| di_smooth_lfhf | 1.26 (700–2800 Hz) | 0.51 (500–2000 Hz) | −0.75 (band follows the point) |
| cardioid_hard / cardioid_aim | 17.15 / 12.36 pass | 17.15 / 12.36 pass | 0 |
| whole_rear (higher is better) | 3.70 | **6.03** | +2.34 |

flat23 falls by 2.5 dB (above the trusted limit, from 10.45 to 0.80 dB, low confidence) because the EQ flattens the listener's direction. Window spread and presence
lift stay put, and the window is still within ±1 dB only up to 713 Hz: EQ changes every direction
alike, so it cannot change how 0–10° differ from 23° at 2–2.85 kHz. The cardioid rows are scored
before the crossover, so they are identical; the whole speaker's rear rejection gains 2.3 dB.
di_smooth_bass is 0.8 dB worse (decision 4).

## Did the round 3 ranking change?

Yes. Screen scores, trusted, in dB. "Sum" is window spread plus presence lift, with round 3's sum in
brackets. Both rank columns count the round 2 reference row (round 3 rank 2, round 3b rank 5). That
row is the same design as the stage base and is not shown.

| 3b rank | Variant | Window | Presence | Sum | DI bass + LF–HF | flat23 | Round 3 rank |
|---|---|---|---|---|---|---|---|
| 1 | bass280-600 | 3.19 | 3.20 | 6.38 (6.71) | 1.19 + 0.83 | 0.85 | 4 |
| 2 | bass260-580 | 3.19 | 3.20 | 6.39 (6.71) | 1.22 + 0.83 | 0.85 | 3 |
| 3 | in5-m310 | 3.17 | 3.19 | 6.36 (6.63) | 1.26 + 0.82 | 1.09 | 5 |
| 4 | stage base (round 2 design) | 3.19 | 3.20 | 6.38 (6.71) | 1.34 + 0.83 | 0.85 | 1 |
| 6 | top60 | 3.40 | 3.29 | 6.69 (6.97) | 1.33 + 0.85 | 0.59 | 8 |
| 7 | top100 | 3.41 | 3.32 | 6.73 (7.07) | 1.37 + 0.84 | 0.90 | 10 |
| 8 | in10-m300 | 3.21 | 3.37 | 6.58 (6.61) | 1.25 + 1.17 | 0.73 | 6 |
| 9 | in10-m290 | 3.22 | 3.42 | 6.64 (6.64) | 1.25 + 1.26 | 0.73 | 7 |
| 10 | in15-m290 | 3.38 | 3.58 | 6.96 (6.76) | 1.26 + 1.24 | 0.73 | 9 |

- **The DI sum orders ranks 1–4.** They share the 6.25–6.50 dB bucket. The raised-bass variants now
  beat the round 2 design by 0.12–0.15 dB of DI. In round 3 the round 2 design won by 0.01 dB.
- **Presence lift improves unevenly.** It improves by 0.25 dB for ranks 1–4, 0.21 dB for top60,
  0.18 dB for top100 and 0.04–0.08 dB for the 10 mm offsets. For in15-m290 it gets 0.11 dB worse.
- **The finals tie.** Both finals fit 170 / 1000 Hz, with filters within 0.05 dB and 0.01 ms, and the
  bass260-580 final equals the full run above. bass280-600 leads only by 0.04 dB of bass DI and
  0.05 dB of whole_rear. The round report recommends it only because it is the first final.

**What this means for the owner's choice:** bass260-580 stands, more firmly than in round 3. It now
leads the round 2 design outright, and raising the woofers another 20 mm buys nothing.

## Decisions for the owner

Each item gives the numbers, then **my recommendation**. The decisions are yours.

### 1. Cap the coax LF boost?

The full-run fit boosts the coax LF by +6.0 dB at 183 Hz (+2.9 dB at 150 Hz). The screen fit boosts
it by +9.5 dB at 150 Hz. The model does not check excursion, so I computed it. Voltages are RMS sine
at the crossover input. SPL is at 2 m and 23°, in the 200–1000 Hz passband.

| Crossover | Coax LF reaches xmax (5.5 mm) | SPL then | One woofer reaches xmax (13.25 mm) | SPL then |
|---|---|---|---|---|
| Fitted, full run (170/1000 Hz) | 87 V at 184 Hz | 121.7 dB | 10.1 V at 20 Hz | **102.9 dB** |
| Fitted, screen (150/1400 Hz) | 51 V at 153 Hz | 117.1 dB | 10.7 V at 20 Hz | 103.6 dB |
| Fixed 170/1400 Hz | 167 V at 201 Hz | 126.9 dB | 10.1 V at 20 Hz | 102.5 dB |

- **The woofers reach xmax first.** They get there about 19 dB before the coax LF with the
  full-run fit, 13.5 dB before with the screen fit and 24.4 dB before with the fixed crossover.
- **The coax LF limit is out of reach.** 87 V at the input is about 100 V at its amplifier, because
  its channel transfer peaks at +1.25 dB. The params assume 40 V.
- **The woofer limit is at 20 Hz, because the bass has no high-pass.** Per 2.83 V, a woofer moves
  3.7 mm at 20 Hz, 3.3 mm at 30 Hz and 2.3 mm at 50 Hz. The sealed F3 is about 44 Hz.

**My recommendation:** keep the fitted boost, and add a subsonic high-pass on the bass at about
30–40 Hz. I did not simulate the high-pass. As an estimate from the excursion per volt: if the
limiting frequency moved from 20 Hz to 30 or 50 Hz, the limit would rise from 10.1 V to 11.3 or
16.5 V (+1.0 to +4.3 dB). This ignores the filter's own slope.

*Method:* lumped-model peak excursion per volt (`cardioid.simulate` for the coax in 16.0 L,
`bass.sealed_excursion` for one woofer in 28.2 L), times each way's electrical transfer over
20–1000 Hz. The SPL is the 2.83 V sum at 23°, scaled. Limits: sine drive, small-signal and linear to
xmax, no thermal, power or amplifier limit, no room gain.

### 2. LF to HF at 1000 Hz

1000 Hz is the bottom of the range; the fixed crossover used 1.4 kHz. It is within an octave of the
compression driver's 750 Hz resonance, and it gives that driver 6–9 dB more drive: −12.4 dB at
1 kHz and −18.0 dB at 750 Hz, against −18.5 and −27.1 dB with the fixed crossover.

Among the fit's candidates (full run), 1000 Hz scores 7.48 dB of window plus presence, 0.51 dB of
LF–HF DI and 1.49 dB of flat23. The best at 1200 Hz scores 7.59, 1.28 and 1.23 dB, and the best at
1400 Hz 8.17, 0.80 and 1.12 dB. So against 1200 Hz, 1000 Hz gains 0.11 dB and 0.77 dB of DI, and
loses 0.26 dB of flat23.

**My recommendation:** check the BMS datasheet's lowest recommended HF crossover, because this model
has no power limit for the compression driver. Keep 1000 Hz if the datasheet allows it; otherwise
use 1200 Hz (−0.11 dB).

### 3. Two filters at the +12 dB limit

The fitter's gain range is −15 to +12 dB, and two filters sit at its top.

**Coax LF, 1893 Hz, Q 0.5.** This is above the coax LF's own 1 kHz low-pass, so the net electrical
level stays low: −1.2 dB at 1 kHz, −7.8 dB at 1.5 kHz and −13.7 dB at 2 kHz, against 0 dB in the
passband. The filter holds up the LR4 slope across the cone's recess notch near 1.9 kHz. That notch
is −11 dB in round 1's computed recess ratio, which rests on estimated cone geometry.

**Coax HF, 12.8 kHz, Q 0.65.** This is far above the 2851 Hz trusted limit, so it rests on the
estimated cone and flare geometry.

**My recommendation:** keep the 1.89 kHz filter, since the cone gets little drive there, and
re-check it once the coax is measured. Leave the 12.8 kHz filter out until there is a measurement,
or set any top-octave lift by ear.

### 4. DI smoothness at the bass–coax crossover got worse

di_smooth_bass rose from 0.55 to 1.35 dB while flat23 fell from 4.01 to 1.49 dB. The ranking did not
trade one against the other: all 36 candidates the fit tried have di_smooth_bass between 1.28 and
1.39 dB. So the loss comes with fitting EQ to 23°, probably the coax LF boost near the bass point.

The 0.25 dB buckets did pick the point. 170/1000 Hz scores 7.480 dB of window plus presence, just
inside the 7.25–7.50 bucket. 190/1000 Hz scores 7.502 dB, which puts it in the next bucket, even
though its DI sum is 0.07 dB lower (1.79 against 1.86). **My recommendation:** accept it. The
alternatives on the grid change it by under 0.1 dB.

### 5. The screen and the full run pick different points

The same design gets 150/1400 Hz on the screen and 170/1000 Hz on the full run. The screen stops at
2715 Hz, so it matches the HF gain over 2.0–2.7 kHz rather than 2–8 kHz. That gives −8.1 dB rather
than −4.8 dB, which shifts the whole LF–HF region. It also starts at 120 Hz rather than 50 Hz, so
the bass EQ and alignment see different data. Neither winner leads by much: 0.14 dB on the screen,
0.02 dB on the full run. **My recommendation:** use the full run's fit (`crossover.txt`), which
sees the whole HF way.

## CAD build

`speaker build` exited 0 in 4 s. This is the first CAD build of the round 2/3 geometry. It wrote 15
STEP, 6 STL and 8 DXF files and `summary.md`, with the guessed parts marked PROVISIONAL. The files
match the design: 1118 × 330 × 350 mm, a 320 mm flare insert, and 12 slots of 70 × 10 mm per side
at a 22 mm pitch, centred 88 mm from the front. The coax chamber is 16.00 L net with 3.68 L of
fillers, the bass is 56.4 L net, and there are no layout warnings.

Two findings, not fixed:
- **The build summary says the cardioid fails.** It gives front/back of 2.0 dB over 150–300 Hz
  ("BELOW the 10 dB target") and 1.1 dB over 300–500 Hz ("BELOW the 6 dB aim"). Its "without
  fillers" block gives 1.6 and 0.8 dB. These are lumped point-source figures, which differ from the
  3D model by up to 28 dB here. The 3D model, which the gate uses, gives 17.15 and 12.36 dB.
- **The directivity report's header mislabels the design.** It calls `params/round-3.toml` "Round 0
  — current design (baseline)". The params path and sha256 under it are correct.

## VituixCAD export

`out/plan3/round-3-design/vituixcad/` (2.9 MB) holds 715 FRD files, 4 ZMA impedance files and
`README.txt` with the import steps. The FRDs are 143 horizontal and vertical polars each for the
coax LF, the coax HF, the bass pair and each woofer alone. **It is SIMULATED, not measured**, and
every file header says so. Levels are for 2 m and 2.83 V, unfiltered. Above 2851 Hz the data are
low confidence, and the bass files stop there.

## Still open

1. **The 2–2.85 kHz window rise.** It is directivity, and neither shape nor EQ has moved it.
2. **Decisions 1–3.**
3. **Driver measurements.** Once bought and measured: refit the crossover, re-size the fillers and
   re-check the recess-notch filter.
4. **M4.** Re-rank the finals on their own scores; it did not matter here, because the finals tie.
