# Round 3b: round 3 re-scored with the fitted crossover

Date: 2026-10-29. Plan 3 replaced the fixed crossover with an active
crossover fitted to each design. The fixed crossover was ideal LR4 at 170 Hz
and 1.4 kHz, with one gain per way. This round reruns round 3 through the
fitted crossover. It also does a full run and a CAD build of the owner's
choice, bass260-580 (`params/round-3.toml`), and exports that design for
VituixCAD. Everything here is simulated, and the drivers are not bought.

## What was run

- **Code commit:** `329132c` (branch `FyrBridd/plan3-crossover-exports`).
- **Round file:** `rounds/round-3.toml`, sha256 `74fa4dd0…`, unchanged from
  round 3.
- **BEM code version:** `d00d81a6…`. Scoring version 2.
- **Commands:**
  1. `uv run speaker directivity --params params/round-3.toml --out out/plan3/round-3-design --trusted-limit 2851`
  2. `uv run speaker round rounds/round-3.toml --out out/rounds/3b`
  3. `uv run speaker build --params params/round-3.toml --out out/plan3/round-3-build`
- **Outputs** (not committed): `out/plan3/round-3-design/` (report,
  `crossover.json`, `crossover.txt`, `vituixcad/`), `out/rounds/3b/report.html`
  and `out/plan3/round-3-build/`.

| Run | Settings | Variants | Errors | Time |
|---|---|---|---|---|
| Full run of bass260-580 | 6 per λ, with the puck | 1 | 0 | 21 s (cached) |
| Round 3b shape stage | screen (5 per λ, to 2715 Hz) | 8 + stage base + round 2 | 0 | 3 min (cached) |
| Round 3b finals | full run | 2 | 0 | 44 min |
| CAD build | – | 1 | 0 | 4 s |

The re-ranked round put bass280-600 into the finals. Its full run had never
been solved, so it was the one new BEM solve (2652 s; about 17 GB of the
30 GB in use at the 4 kHz band). Everything else came from the cache.

## The fitted crossover of bass260-580

Points: **bass to coax 170 Hz, coax LF to HF 1000 Hz**. The search ranges
were 150–250 Hz and 1.0–2.0 kHz. The acoustic target is LR4 at each point,
with the 23° response flat.

| Way | Gain | Delay | Polarity | PEQ (f0, gain, Q) |
|---|---|---|---|---|
| bass | +11.5 dB | 1.25 ms | inverted | 80 Hz +1.2 dB Q1.8; 213 Hz −1.8 dB Q1.3 |
| coax LF | 0 dB | 0 | normal | 183 Hz +6.0 dB Q2.5; 585 Hz −5.3 dB Q0.9; 1394 Hz −2.7 dB Q3.4; 1893 Hz +12.0 dB Q0.5 |
| coax HF | −4.8 dB | 0.01 ms | normal | 1885 Hz −4.7 dB Q1.0; 3997 Hz −6.6 dB Q2.9; 5620 Hz +1.9 dB Q5.4; 11831 Hz +2.0 dB Q3.2; 12833 Hz +12.0 dB Q0.6 |

- The delays are well inside their search grids (±3 ms and ±1 ms).
- At 23°, the sum stays within −1.5 / +1.6 dB of its 200–1000 Hz level over
  150–2851 Hz. That level is 91.9 dB at 2 m for 2.83 V. The sum has no
  notches. The fixed crossover gave −3.9 / +4.1 dB.
- `crossover.txt` is the DSP setup sheet: these values, plus biquads at 48
  and 96 kHz.

**Passive-network warnings:**
1. The coax LF needs 10.9 dB of summed EQ boost at 2 kHz. A passive network
   can supply 6 dB.
2. The coax HF needs 13.7 dB of summed EQ boost at 12.1 kHz.
3. The 1000 Hz point is within an octave of the 750 Hz HF resonance (the
   limit is 1500 Hz).
4. The coax HF filters at 4.0, 5.6, 11.8 and 12.8 kHz are above the trusted
   limit (2851 Hz), so they rest on low-confidence data.

## Scores against the fixed-crossover final

bass260-580, full run, in dB, trusted unless marked:

| Score | Fixed (round 3 final) | Fitted | Change |
|---|---|---|---|
| flat23 (150–2851) | 4.01 | **1.49** | −2.52 |
| flat23 (above, low) | 10.45 | 0.80 | −9.64 |
| window_spread (150–2851) | 3.48 warn | 3.36 warn | −0.13 |
| presence_lift (2000–2851) | 4.10 warn | 4.12 warn | +0.02 |
| presence_lift (2851–5000, low) | 6.08 | 5.54 | −0.54 |
| di_smooth_bass | 0.55 | **1.35** | +0.80 (worse) |
| di_smooth_lfhf | 1.26 (700–2800 Hz) | 0.51 (500–2000 Hz) | −0.75 (band follows the point) |
| cardioid_hard (150–300) | 17.15 pass | 17.15 pass | 0 |
| cardioid_aim (300–500) | 12.36 pass | 12.36 pass | 0 |
| whole_rear (higher is better) | 3.70 | **6.03** | +2.34 |

- **flat23 falls by 2.5 dB.** That was expected: the EQ flattens the 23°
  response, which is the listener's direction.
- **Window spread and presence lift do not move.** The 0–10° directions
  still rise over 23° at 2–2.85 kHz. The window is within ±1 dB only up to
  713 Hz, as before. EQ changes every direction alike, so it cannot change
  how the directions differ from each other.
- **The cardioid rows are identical**, because they score the coax section
  before the crossover. After the crossover, the whole speaker's rear
  rejection gains 2.3 dB.
- **DI smoothness at the bass–coax crossover is 0.8 dB worse** (decision 4).

## Did the round 3 ranking change?

Yes. Screen scores, trusted, in dB. The sum is window spread plus presence
lift, with round 3's sum in brackets.

| 3b rank | Variant | Window | Presence | Sum | DI bass + LF–HF | flat23 | Round 3 rank |
|---|---|---|---|---|---|---|---|
| 1 | bass280-600 | 3.19 | 3.20 | 6.38 (6.71) | 1.19 + 0.83 | 0.85 | 3 |
| 2 | bass260-580 | 3.19 | 3.20 | 6.39 (6.71) | 1.22 + 0.83 | 0.85 | 2 |
| 3 | in5-m310 | 3.17 | 3.19 | 6.36 (6.63) | 1.26 + 0.82 | 1.09 | 4 |
| 4 | stage base (round 2 design) | 3.19 | 3.20 | 6.38 (6.71) | 1.34 + 0.83 | 0.85 | 1 |
| 6 | top60 | 3.40 | 3.29 | 6.69 (6.97) | 1.33 + 0.85 | 0.59 | 7 |
| 7 | top100 | 3.41 | 3.32 | 6.73 (7.07) | 1.37 + 0.84 | 0.90 | 9 |
| 8 | in10-m300 | 3.21 | 3.37 | 6.58 (6.61) | 1.25 + 1.17 | 0.73 | 5 |
| 9 | in10-m290 | 3.22 | 3.42 | 6.64 (6.64) | 1.25 + 1.26 | 0.73 | 6 |
| 10 | in15-m290 | 3.38 | 3.58 | 6.96 (6.76) | 1.26 + 1.24 | 0.73 | 8 |

Rank 5 is the round 2 reference, which is the same design as the stage base.

- Ranks 1–4 share the 6.25–6.50 dB bucket, so the DI sum orders them. The
  raised-bass variants now beat the round 2 design by 0.12–0.15 dB of DI.
  In round 3 the round 2 design won by 0.01 dB.
- Most variants improve presence lift by 0.2–0.25 dB on the screen.
- The coax offset variants gain least: 0.04–0.08 dB, and at 15 mm presence
  lift gets 0.11 dB worse. Offsetting the coax is not worth a handed pair.

**Finals.** Both finals fit the same crossover, 170 / 1000 Hz. Every gain,
delay and filter agrees within 0.05 dB and 0.01 ms. The bass260-580 final
reproduces the full run above exactly. bass280-600 leads by 0.04 dB of
bass DI smoothness (1.31 against 1.35 dB) and by 0.05 dB of whole-speaker
rear rejection (6.08 against 6.03 dB). Every other trusted and
low-confidence row agrees to 0.01 dB. The round report recommends
bass280-600 because it is the first final.

**What this means for the owner's choice:** bass260-580 stands, and the
case for it is stronger than in round 3. With the fitted crossover it beats
the round 2 design on the screen by 0.12 dB of DI, where round 3 had it
behind by a 0.01 dB tie-break. bass280-600 is a tie (0.04 dB), so there is
no reason to raise the woofers another 20 mm or change the drawings.

## Decisions for the owner

Each item gives the numbers, then **my recommendation**. The decisions are
yours.

### 1. Cap the coax LF boost?

The fit boosts the coax LF just above the bass–coax point. In the full run
the boost is +6.0 dB at 183 Hz, which is +2.9 dB at 150 Hz. The screen fit
of the same design boosts +9.5 dB at 150 Hz. The model does not check
excursion, so I computed it separately:

| Crossover | Coax LF reaches xmax (5.5 mm) | System SPL then | One woofer reaches xmax (13.25 mm) | System SPL then |
|---|---|---|---|---|
| Fitted, full run (170/1000 Hz) | 87 V at 184 Hz | 121.7 dB | 10.1 V at 20 Hz | **102.9 dB** |
| Fitted, screen (150/1400 Hz) | 51 V at 153 Hz | 117.1 dB | 10.7 V at 20 Hz | 103.6 dB |
| Fixed 170/1400 Hz | 167 V at 201 Hz | 126.9 dB | 10.1 V at 20 Hz | 102.5 dB |

Voltages are RMS sine at the crossover input. SPL is at 2 m and 23°, in the
200–1000 Hz passband.

- **The woofers reach xmax 13–14 dB before the coax LF does**, with either
  fit. The boost costs the coax 5–10 dB of headroom against the fixed
  crossover, but the coax LF is still far from its limit. At its full-fit
  limit, the coax amplifier would have to deliver about 100 V; the params
  assume 40 V.
- The woofers' limit is at 20 Hz, because the bass way has no high-pass.
  Per 2.83 V of input, a woofer moves 3.7 mm at 20 Hz, 3.3 mm at 30 Hz and
  2.3 mm at 50 Hz. The sealed bass's F3 is about 44 Hz, so little output
  comes from that range.
- **My recommendation: keep the fitted boost, and add a subsonic high-pass
  on the bass at about 30–40 Hz.** By the figures above, that should raise
  the excursion limit by roughly 1–4 dB. I did not simulate the
  high-pass. Raising the bass–coax point to unload the coax would not help,
  because the coax is not the limit.

**Method and limits.**
- The peak cone excursion per volt comes from the design's lumped models.
  The coax uses `speaker.cardioid.simulate` in its 16.0 L vented chamber.
  Each woofer uses `speaker.bass.sealed_excursion` in half of 56.4 L. The
  woofers are in parallel, so each gets the full bass channel voltage.
- I multiplied that by each way's electrical crossover transfer over
  20–1000 Hz. The worst frequency gives the voltage at xmax.
- The SPL is the 2.83 V crossover sum at 23°, scaled to that voltage.
- Limits: sine drive; a small-signal lumped model, linear up to xmax; no
  thermal, power or amplifier limit; no room gain. The woofer column is
  the passband SPL at the moment the woofer hits xmax at 20 Hz, not the SPL
  at 20 Hz.

### 2. LF to HF at 1000 Hz

The fit chose 1000 Hz, the bottom of its range, where the fixed crossover
used 1.4 kHz. The costs:
- 1000 Hz is within an octave of the compression driver's 750 Hz
  resonance (passive warning 3).
- The compression driver gets more drive. Its electrical level is −12.4 dB
  at 1 kHz and −18.0 dB at 750 Hz. The fixed crossover gave −18.5 and
  −27.1 dB, so this is 6–9 dB more near the resonance.

The fit's own candidates, in dB (full run, best bass point for each
LF–HF point):

| LF–HF | Window + presence | DI bass + LF–HF | flat23 |
|---|---|---|---|
| **1000 Hz** (bass 170 Hz, chosen) | 7.48 | 1.35 + 0.51 | 1.49 |
| 1200 Hz (bass 150 Hz) | 7.59 | 1.29 + 1.28 | 1.23 |
| 1400 Hz (bass 150 Hz) | 8.17 | 1.28 + 0.80 | 1.12 |
| 1600 Hz (bass 150 Hz) | 8.42 | 1.28 + 1.12 | 0.99 |

Against 1200 Hz, 1000 Hz gains 0.11 dB of window plus presence and 0.7 dB of
LF–HF DI smoothness, and loses 0.26 dB of flat23. Against 1400 Hz it gains
0.7 dB of window plus presence.

**My recommendation:** check the BMS datasheet's lowest recommended
crossover for the HF section before using 1000 Hz. This model has no power
or excursion limit for the compression driver. If the datasheet allows
1 kHz, keep it. If not, use 1200 Hz, which costs 0.11 dB of window plus
presence. 1200 Hz is still within an octave of the resonance, so check it
against the datasheet too.

### 3. Two filters at the +12 dB limit

The fitter's gain range is −15 to +12 dB, and two filters sit at its top.

- **Coax LF, +12.0 dB at 1893 Hz, Q 0.5.** This is above the coax LF's own
  1 kHz low-pass, so the net electrical level stays low: −1.2 dB at 1 kHz,
  −7.8 dB at 1.5 kHz, −13.7 dB at 2 kHz and −21.7 dB at 2.5 kHz, against
  0 dB in the passband. The filter holds up the LR4 slope across the
  cone's recess notch near 1.9 kHz. That notch is −11 dB in the recess
  ratio round 1 computed, and it rests on the estimated cone geometry.
- **Coax HF, +12.0 dB at 12.8 kHz, Q 0.65.** This is far above the 2851 Hz
  trusted limit. There the response rests on the estimated cone and flare
  geometry.

**My recommendation:** keep the 1.89 kHz filter, since the cone gets
little drive there. Re-check it once the coax is measured, because it is
fitted to a modelled notch. Leave the 12.8 kHz filter out of the DSP, or
set any top-octave lift by ear, until there is a measurement.

### 4. DI smoothness at the bass–coax crossover got worse

di_smooth_bass rose from 0.55 to 1.35 dB while flat23 fell from 4.01 to
1.49 dB.
- **The ranking did not trade the two.** All 36 candidates the fit tried
  have di_smooth_bass between 1.28 and 1.39 dB, at every bass–coax point.
  So the loss comes with fitting EQ to the 23° response, probably the coax
  LF boost near the bass point. The choice of point did not cause it.
- **The buckets act only at the margin.** The search sorts candidates into
  0.25 dB buckets of window plus presence, and the DI sum decides only
  within a bucket. 170/1000 Hz scores 7.480 dB, just inside the 7.25–7.50
  bucket. 190/1000 Hz scores 7.502 dB, so it falls into the next bucket,
  although its DI sum is 0.07 dB lower (1.79 against 1.86). A 0.02 dB edge
  in window plus presence beat a 0.07 dB edge in DI.

DI smoothness is about how evenly the radiated power changes through the
crossover. That shapes the room sound more than the direct sound at the
seat.

**My recommendation:** accept it. The alternatives within the fit's grid
change it by under 0.1 dB.

### 5. The screen and the full run pick different points

The same design gets 150/1400 Hz on the screen and 170/1000 Hz on the full
run. The grids differ:
- **The screen stops at 2715 Hz.** It matches the HF gain only over
  2.0–2.7 kHz, where the full run uses 2–8 kHz. The HF gain comes out at
  −8.1 dB on the screen and −4.8 dB on the full run, which moves the whole
  LF–HF region.
- **The screen starts at 120 Hz, the full run at 50 Hz.** So the bass EQ
  and alignment also see different data.

Neither winner is a strong preference. The screen's winner leads the
runner-up by 0.14 dB of window plus presence (6.39 against 6.52 dB). The
full run's winner leads by 0.02 dB.

**My recommendation:** use the full run's fit (`crossover.txt`), which sees
the whole HF way. Treat the screen's points as a ranking aid only. The
ranking does not hang on them: the screen's top four are within 0.03 dB on
window plus presence.

## CAD build

`speaker build` exit 0, 4 s. This is the first CAD build of the round 2/3
geometry. It wrote 15 STEP, 6 STL and 8 DXF files and `summary.md`, with
the guessed parts marked PROVISIONAL. Checked in the files:
- the cabinet is 1118 × 330 × 350 mm;
- the flare insert is 320 mm across;
- each side has 12 slots of 70 × 10 mm at a 22 mm pitch, centred 88 mm
  from the front;
- the coax chamber is 16.00 L net with 3.68 L of fillers, and the bass is
  56.4 L net;
- there are no layout warnings.

Two findings, not fixed:
- **The build summary reports the cardioid as failing.** Its line reads
  "Minimum front/back over 150–300 Hz: 2.0 dB (BELOW the 10 dB target)".
  That figure is the lumped point-source model, which the directivity
  report says differs from the 3D model by up to 28 dB here. The 3D
  model, which is what the gate uses, gives 17.15 dB.
- **The directivity report's header mislabels the design.** It calls
  `params/round-3.toml` "Round 0 — current design (baseline)". The params
  path and sha256 under the header are correct.

## VituixCAD export

`out/plan3/round-3-design/vituixcad/` (2.9 MB). It has 715 FRD files: 143
polar files, horizontal and vertical, for each of five sets (coax LF, coax
HF, the bass pair, and each woofer alone). It also has 4 ZMA impedance
files and a `README.txt` with the import steps.

**It is SIMULATED, not measured**, and every file header says so. Levels
are for 2 m and 2.83 V, unfiltered. Above 2851 Hz the data are low
confidence, and the bass files stop there. No row is non-finite.

## Still open

1. **The 2–2.85 kHz window rise** (presence lift 4.1 dB). It is a
   directivity problem, so the EQ could not touch it. Neither the shape
   rounds nor the crossover has moved it.
2. **Decisions 1–3**: the bass high-pass, the HF point and the two +12 dB
   filters.
3. **The drivers are not bought.** Once they are measured: refit the
   crossover, re-size the fillers, and re-check the recess-notch filter.
4. **Re-rank the finals on their own scores (M4)**, carried over from round
   3. Here it did not matter, as the finals tie.
