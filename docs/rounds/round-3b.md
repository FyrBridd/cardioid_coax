# Round 3b: round 3 re-scored with the fitted crossover

Date: 2026-10-29. Plan 3 replaced the fixed crossover (ideal LR4 at 170 Hz and 1.4 kHz, one gain per
way) with an active crossover fitted to each design. This round reruns round 3 through it. It also
makes a full run, a CAD build and a VituixCAD export of the owner's choice, bass260-580
(`params/round-3.toml`). Everything here is simulated; the drivers are not bought.

## What was run

- **Versions:** code commit `2795be3` (branch `FyrBridd/plan3-crossover-exports`), BEM code version
  `d00d81a6…`, scoring version 2. The round file `rounds/round-3.toml` (sha256 `74fa4dd0…`) is
  unchanged. Outputs are not committed.
- **This is the second version of this write-up.** The first used a crossover search that read data
  above the trusted limit (below, decision 5). The search now uses the trusted band only, and the
  full run and the finals were rerun with it.
- **Full run:** `uv run speaker directivity --params params/round-3.toml --out out/plan3/round-3-design --trusted-limit 2851`
  took 16 s, all cached.
- **Round:** `uv run speaker round rounds/round-3.toml --out out/rounds/3b` took 31 s. The shape stage
  was reused (10 variants, 0 errors); the finals took 14 s each, cached. Re-scoring two screen variants
  with the new code gave identical scores, so the screen is unchanged.
- **CAD:** `uv run speaker build --params params/round-3.toml --out out/plan3/round-3-build` took 4 s.
- **The one new BEM solve** of round 3b was the bass280-600 final, in the first version (2652 s).

## The fitted crossover of bass260-580

Points: **bass to coax 230 Hz, coax LF to HF 1000 Hz** (search ranges 150–250 Hz and 1.0–2.0 kHz).
The acoustic target is LR4 at each point, flat at 23°. Every choice the fit makes uses the data up to
the trusted limit (2851 Hz). The filters marked *low* come after the search, from the data above it.
They move the response below 2851 Hz by at most 0.25 dB.

| Way | Gain | Delay | Polarity | PEQ (f0, gain, Q) |
|---|---|---|---|---|
| bass | +11.5 dB | 1.27 ms | inverted | 82 Hz +1.1 dB Q2.3; 212 Hz −1.8 dB Q1.6 |
| coax LF | 0 dB | 0.02 ms | normal | 236 Hz +2.6 dB Q6.0; 594 Hz −5.2 dB Q1.0; 1392 Hz −2.8 dB Q3.1; 1893 Hz +12.0 dB Q0.5 (limit) |
| coax HF | −8.0 dB | 0 | normal | 1101 Hz +1.8 dB Q0.5; 1748 Hz −2.6 dB Q2.0; *low:* 4308 Hz −9.4 dB Q0.8; 6401 Hz +12.0 dB Q0.63 (limit); 13188 Hz +12.0 dB Q1.1 (limit) |

The previous version had 170 / 1000 Hz and a coax HF gain of −4.8 dB, matched over 2–8 kHz rather
than 2.0–2.85 kHz. The delays are well inside their grids (±3 ms and ±1 ms). Over 150–2851 Hz the
sum at 23° stays within −0.8 / +1.1 dB of its 200–1000 Hz level (91.75 dB at 2 m for 2.83 V), with
no notches. The fixed crossover gave −3.9 / +4.1 dB, the previous fit −1.5 / +1.6.

`crossover.txt` is the DSP setup sheet. Its header gives the design file and sha256, the trusted
limit, the largest gain plus PEQ boost of any way (+12.5 dB, the bass at 81 Hz, for setting gains
that do not clip) and a note that miniDSP-style tools want −a1 and −a2.

**48 or 96 kHz.** The biquads are RBJ sections. How far each cascade misses its analog way (largest
|dB| where the way is within 20 dB of its peak):

| Way | 48 kHz: trusted / 2–5 kHz / above 5 kHz | 96 kHz: trusted / 2–5 kHz / above 5 kHz |
|---|---|---|
| bass | 0.00 / – / – | 0.00 / – / – |
| coax LF | 0.26 / 0.26 / – | 0.06 / 0.06 / – |
| coax HF | 0.50 / **1.03** / **7.01** | 0.13 / 0.26 / 1.61 |

At 48 kHz the bilinear transform squeezes the broad HF filters at 6.4 and 13.2 kHz, so the sheet
says **use 96 kHz** (decision 3). Without the three *low* HF filters, the 48 kHz coax HF error falls
to 0.03 dB; the coax LF's 0.26 dB (its 1.9 kHz boost) stays.

**Passive-network warnings:**
1. EQ boost: the coax LF needs 11.0 dB at 2 kHz and the coax HF 16.4 dB at 12.7 kHz (passive: 6 dB).
2. The 1000 Hz point is within an octave of the 750 Hz HF resonance (the limit is 1500 Hz).
3. The coax HF filters at 4.3, 6.4 and 13.2 kHz are above the trusted limit (2851 Hz).

## Scores against the fixed-crossover final

bass260-580, full run, in dB, trusted unless marked:

| Score | Fixed (round 3 final) | Previous fit (170/1000) | Fitted (230/1000) |
|---|---|---|---|
| flat23 (150–2851 Hz) | 4.01 | 1.49 | **0.92** |
| flat23 (2851–16000 Hz, low confidence) | 10.45 | 0.80 | 2.45 |
| window_spread (150–2851 Hz) | 3.48 warn | 3.36 warn | 3.45 warn |
| presence_lift (2000–2851 Hz) | 4.10 warn | 4.12 warn | **3.64 warn** |
| presence_lift (2851–5000 Hz, low confidence) | 6.08 | 5.54 | 4.41 |
| di_smooth_bass | 0.55 (85–340 Hz) | 1.35 (85–340 Hz) | **1.33** (115–460 Hz) |
| di_smooth_lfhf | 1.26 (700–2800 Hz) | 0.51 (500–2000 Hz) | 0.39 (500–2000 Hz) |
| cardioid_hard / cardioid_aim | 17.15 / 12.36 pass | same | same |
| whole_rear (higher is better) | 3.70 | 6.03 | 4.16 |

Against the fixed crossover, flat23 falls by 3.1 dB and presence lift by 0.46 dB (0.48 dB below the
previous fit, whose HF gain was 3.2 dB higher). Window spread does not move; the window stays within
±1 dB only up to 713 Hz, because EQ changes every direction alike. The cardioid rows are scored
before the crossover. Whole-speaker rear rejection is 0.46 dB better than fixed but 1.9 dB below the
previous fit, as the bass point moved to 230 Hz. di_smooth_bass is 0.8 dB worse (decision 4).

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
- **The finals tie.** Both fit 230 / 1000 Hz with filters within 0.05 dB and 0.01 ms; bass260-580's
  equals the full run above. bass280-600 leads by 0.07 dB of bass DI (1.26 against 1.33) and 0.04 dB
  of whole_rear, other trusted scores within 0.01 dB. The report recommends it only as the first final.
  The screen ranking is unchanged by the trusted-only search: the screen grid stops at the limit.

**What this means for the owner's choice:** bass260-580 stands, more firmly than in round 3. It now
leads the round 2 design outright, and raising the woofers another 20 mm buys nothing.

## Decisions for the owner

Each item gives the numbers, then **my recommendation**. The decisions are yours.

### 1. Cap the coax LF boost?

The full-run fit now boosts the coax LF by only +2.6 dB, at 236 Hz (Q 6). The previous fit boosted it
by +6.0 dB at 183 Hz, and the screen fit by +9.5 dB at 150 Hz. The model does not check excursion,
so I computed it. Voltages are RMS sine at the crossover input. SPL is at 2 m and 23°, in the
200–1000 Hz passband.

| Crossover | Coax LF reaches xmax (5.5 mm) | SPL then | One woofer reaches xmax (13.25 mm) | SPL then |
|---|---|---|---|---|
| Fitted, full run (230/1000 Hz) | 198 V at 237 Hz | 128.6 dB | 10.1 V at 20 Hz | **102.8 dB** |
| Fitted, screen (150/1400 Hz) | 51 V at 153 Hz | 117.1 dB | 10.7 V at 20 Hz | 103.6 dB |
| Fixed 170/1400 Hz | 167 V at 201 Hz | 126.9 dB | 10.1 V at 20 Hz | 102.5 dB |

- **The woofers reach xmax first,** about 26 dB before the coax LF with the full-run fit. A cap would
  change nothing: a +3 dB cap leaves the fit as it is.
- **The woofer limit is at 20 Hz, because the bass has no high-pass.** Per 2.83 V, a woofer moves
  3.7 mm at 20 Hz, 3.3 mm at 30 Hz and 2.3 mm at 50 Hz. The sealed F3 is about 44 Hz.

**My recommendation:** no cap. Add a subsonic high-pass on the bass at about 30–40 Hz. I did not
simulate the high-pass. As an estimate from the excursion per volt: if the limiting frequency moved
from 20 Hz to 30 or 50 Hz, the limit would rise from 10.1 V to 11.3 or 16.6 V (+1.0 to +4.3 dB).
This ignores the filter's own slope.

*Method:* lumped-model peak excursion per volt (`cardioid.simulate` for the coax in 16.0 L,
`bass.sealed_excursion` for one woofer in 28.2 L), times each way's electrical transfer over
20–1000 Hz. The SPL is the 2.83 V sum at 23°, scaled. Limits: sine drive, small-signal and linear to
xmax, no thermal, power or amplifier limit, no room gain.

### 2. LF to HF at 1000 Hz

1000 Hz is the bottom of the range; the fixed crossover used 1.4 kHz. It is within an octave of the
compression driver's 750 Hz resonance, and it gives that driver 6–8 dB more drive: −12.5 dB at
1 kHz and −18.9 dB at 750 Hz, against −18.5 and −27.1 dB with the fixed crossover.

Among the fit's 36 candidates (scored on the trusted band), 230/1000 Hz scores 7.16 dB of window plus
presence, 0.44 dB of LF–HF DI and 1.17 dB of flat23. The best at 1200 Hz (150/1200) scores 7.13, 1.25
and 0.89 dB, and the best at 1400 Hz (150/1400) 7.63, 0.69 and 0.96 dB. So against 1200 Hz, 1000 Hz
loses 0.03 dB of window plus presence and 0.28 dB of flat23, and gains 0.81 dB of LF–HF DI.

**My recommendation:** check the BMS datasheet's lowest recommended HF crossover, because this model
has no power limit for the compression driver. Keep 1000 Hz if the datasheet allows it; otherwise
use 1200 Hz; the fit then picks 150 Hz for the bass point, and the window plus presence is no worse.

### 3. Three filters at the +12 dB limit, and the sample rate

The fitter's gain range is −15 to +12 dB. Three filters sit at its top; the sheet and report mark them.

**Coax LF, 1893 Hz, Q 0.5.** This is above the coax LF's own 1 kHz low-pass, so the net electrical
level stays low: −1.2 dB at 1 kHz, −7.8 dB at 1.5 kHz and −13.6 dB at 2 kHz, against 0 dB in the
passband. The filter holds up the LR4 slope across the cone's recess notch near 1.9 kHz. That notch
is −11 dB in round 1's computed recess ratio, which rests on estimated cone geometry.

**Coax HF, 6.4 kHz (Q 0.63) and 13.2 kHz (Q 1.1),** with a −9.4 dB filter at 4.3 kHz. All three are
above the 2851 Hz trusted limit, so they rest on the estimated cone and flare geometry. They are
also what makes 48 kHz miss by 1 dB.

**My recommendation:** keep the 1.89 kHz filter, since the cone gets little drive there, and
re-check it once the coax is measured. Leave the three HF filters above 2.85 kHz out until there is
a measurement, or set any top-octave lift by ear. Run the DSP at 96 kHz if it can. On a 48 kHz-only
DSP, leaving out those three filters brings every way within 0.26 dB of its analog filter.

### 4. DI smoothness at the bass–coax crossover got worse

di_smooth_bass rose from 0.55 to 1.33 dB while flat23 fell from 4.01 to 0.92 dB. The ranking did not
trade one against the other: all 36 candidates have di_smooth_bass between 1.28 and 1.39 dB. So the
loss comes with fitting EQ to 23°.

230/1000 Hz leads its 0.25 dB bucket of window plus presence (7.00–7.25 dB) on the DI sum: 1.77 dB,
against 1.81 for 210/1000. 150–190/1000 Hz fall two buckets (8.30–8.36 dB). There the HF alignment
puts the HF 0.46 ms ahead of the LF, which raises the window spread to 4.29 dB. **My
recommendation:** accept 230 Hz. Its neighbours at 210 and 250 Hz score within 0.06 dB.

### 5. The screen and the full run still pick different points

The previous version recommended the full-run fit because it "sees the whole HF way". That was a
mistake: its HF gain was matched over 2–8 kHz and its EQ fitted up to 19 kHz, so the point search
rested on low-confidence data, against spec §4. Both fits now use only the data up to 2851 Hz, and
their HF gains agree (−8.0 dB full run, −8.1 dB screen). The points still differ, 150/1400 Hz on the
screen and 230/1000 Hz on the full run: the screen starts at 120 Hz rather than 50 Hz and uses a
coarser mesh (5 rather than 6 per wavelength). Neither winner leads by much: 0.14 dB on the screen,
0.05 dB of DI on the full run. **My recommendation:** use the full run's fit (`crossover.txt`): its
data are finer, and its search now uses the same trusted band as the screen's.

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
low confidence, and the bass files stop there. The README says to set VituixCAD's listening
distance to 2 m or your seat, and why.

## Still open

1. **The 2–2.85 kHz window rise.** It is directivity, and neither shape nor EQ has moved it.
2. **Decisions 1–3.**
3. **Driver measurements.** Once bought and measured: refit the crossover, re-size the fillers and
   re-check the recess-notch filter.
4. **M4.** Re-rank the finals on their own scores; it did not matter here, because the finals tie.
