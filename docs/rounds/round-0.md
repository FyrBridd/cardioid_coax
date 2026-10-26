# Round 0: directivity baseline

Date: 2026-10-26. This is the current design, scored as it stands. Nothing
was optimised.

## What was run

- **Params commit:** `21b7d79` (branch `FyrBridd/Kickoff`).
- **Params file:** `params/default.toml`, sha256
  `f877a7c7c80c008681db9059f1a5f6a65ba9631fb23489afc5800a8cdc86b887`.
- **BEM code version** (the cache key):
  `f6c88fea0e27ca2004464df34e727462f76f43cee9b2ef99b050ca46d34be3ea`.
- **Commands:**
  1. `uv run speaker directivity --mesh-check --out out/directivity`
  2. `uv run speaker directivity --out out/directivity`
- **Outputs** (not committed; `out/` is a build folder):
  `out/directivity/report.html`, `scores.json`, `mesh_check.json` and
  `plots/*.png`.
- **Machine:** 16 cores, 30 GB RAM, no GPU. The BEM memory budget is 22 GB.

## Trusted 3D limit: 2.85 kHz

**How it was measured.** The mesh check solves the cabinet twice: once with
the full run's mesh (6 elements per wavelength) and once with a finer mesh
(8 per wavelength). The finer mesh also shrinks the three size caps by 6/8:

| Mesh size cap | Full run | Check mesh |
|---|---|---|
| Whole cabinet (h_max) | 50 mm | 37.5 mm |
| Near the coax cavity (h_coax) | 10 mm | 7.5 mm |
| Near the vents (h_vent) | 20 mm | 15 mm |

It then compares these curves between the two meshes, per 1/3 octave:

- the 23° curve, listening window, early reflections and sound power;
- H 0–30° in 5° steps, and H ±40–90° in 10° steps, after the crossover;
- the coax section's rear rejection, L(0°) − L(θ) for θ = ±120–180°, over
  150–500 Hz.

A band is trusted if every curve agrees within 0.5 dB, and every band below
it does too.

**Result.** Every band that was checked agreed. The worst difference was
0.08 dB (the 160 Hz band). No band came near the 0.5 dB limit.

**Why it stops at 2.85 kHz.** Memory, not accuracy. The check solved bands up
to 3000 Hz. The highest frequency on the grid in that band is 2851 Hz. The
next band (to 4000 Hz) needs a check mesh of about 24.0 GB. The budget is
22.0 GB, so the check stopped there. So 2.85 kHz is the highest frequency
the check could verify on this machine. The true limit of the full-run mesh
may be higher. We cannot show it here.

## Scores

All 15 rows of `scores.json`. "Low" rows are above the trusted limit. They
are shown, but they do not rank designs.

"Puck share" is the share of the band's frequency points that use any of the
puck model (the coax on a round cylinder, with no cabinet), and the mean puck
weight over the band. See observation 3.

| Key | Score | Band (Hz) | Confidence | Value | Unit | Status | Worst at | Puck share |
|---|---|---|---|---|---|---|---|---|
| flat23 | 23° H flatness | 150–2851 | trusted | 4.10 | dB | info | – (RMS 2.87 dB) | 24% / 0.13 |
| flat23 | 23° H flatness | 2851–16000 | low | 9.57 | dB | info | – (RMS 6.18 dB) | 100% / 1.00 |
| window_spread | Window spread vs 23° | 150–2851 | trusted | 4.70 | dB | warn | within ±1 dB only from 150 to 898 Hz | 24% / 0.13 |
| window_spread | Window spread vs 23° | 2851–16000 | low | 15.37 | dB | warn | – | 100% / 1.00 |
| presence_lift | Presence lift 2–5 kHz | 2000–2851 | trusted | 5.67 | dB | warn | – | 100% / 0.72 |
| presence_lift | Presence lift 2–5 kHz | 2851–5000 | low | 7.15 | dB | warn | – | 100% / 1.00 |
| di_smooth_bass | DI smoothness at bass–coax crossover | 85–340 | trusted | 0.54 | dB | info | SPDI 0.54, ERDI 0.12 | 0% |
| di_smooth_lfhf | DI smoothness at LF–HF crossover | 700–2800 | trusted | 1.42 | dB | info | SPDI 1.42, ERDI 0.72 | 48% / 0.24 |
| rolloff_monotonic | Off-axis rolloff monotonic | 150–2851 | trusted | 0 | count | pass | – | 24% / 0.13 |
| rolloff_monotonic | Off-axis rolloff monotonic | 2851–16000 | low | 2 | count | warn | 12.7 and 16 kHz, both sides | 100% / 1.00 |
| cardioid_hard | Cardioid rejection 150–300 Hz | 150–300 | trusted | 5.28 | dB | **fail** | 300 Hz, 180° | 0% |
| cardioid_hard_far | Cardioid rejection 150–300 Hz at 50 m | 150–300 | trusted | 5.57 | dB | info | 300 Hz | 0% |
| cardioid_aim | Cardioid rejection 300–500 Hz | 300–500 | trusted | 5.13 | dB | warn | 336 Hz, 180° | 0% |
| rear_hemisphere | Rear hemisphere rejection | 150–500 | trusted | 5.13 | dB | info | 336 Hz, 180° | 0% |
| whole_rear | Whole-system rear rejection | 150–300 | trusted | 2.71 | dB | info | 150 Hz, 180° | 0% |

Vent resistance spread (the felt as built is unknown, ±30% around
R = 190 Pa·s/m):

| | R × 0.7 | R × 1.0 | R × 1.3 |
|---|---|---|---|
| cardioid_hard (target ≥ 10 dB) | 3.5 dB | 5.3 dB | 6.2 dB |
| cardioid_aim (target ≥ 6 dB) | 3.7 dB | 5.1 dB | 5.7 dB |

Crossover gains (fixed LR4 at 170 Hz and 1400 Hz, no EQ, no delays; each way
matched at 23° to the coax LF way):

| Way | Gain | Matching band |
|---|---|---|
| bass | +14.19 dB | 80–150 Hz |
| coax LF | 0.00 dB | 300–1000 Hz |
| coax HF | −2.32 dB | 2000–8000 Hz |

## The three things that matter most

### 1. The cardioid fails its hard gate

The coax section rejects only 5.3 dB to the rear over 150–300 Hz. The target
is at least 10 dB. The worst point is 300 Hz, straight behind (180°).

It is not a near-field effect of the 2 m distance. At 50 m the same number
is 5.6 dB.

The vent felt probably will not fix it on its own:

- Across the ±30% felt spread the score moves from 3.5 to 6.2 dB.
- From the Task 7 investigation (2026-10-26), not reproduced in this run's
  outputs: sweeps of the vent resistance with the 3D fields found that
  R ≈ 400–500 Pa·s/m gives about 13 dB at 150 Hz, but only about 8.6 dB at
  250 Hz. (The current R is 190 Pa·s/m.)

So plan 2 should probably look at the vent position, area and shape, not
just the felt.

**Why the simple model was so far off.** The lumped (free-field) model
differs from the 3D result by up to 19.4 dB over 150–500 Hz (report). The
following come from the Task 7 investigation (2026-10-26), not reproduced in
this run's outputs: the lumped model predicted 22–24 dB in the hard band, and
there are three likely causes of the gap:

- The cabinet shadows and diffracts the vent sound. The vents seem to act as
  if they were about 0.43 m from the cone, not 0.234 m.
- In 3D the cone radiates from the mouth of its cavity, not from the baffle
  plane the lumped model assumes.
- The score is taken at 2 m, not in the far field.

**What the vents do now.** I split the coax section into its two parts using
the cached 3D fields (the cone's term and the vents' term of the same
solution, at 2 m, before the crossover). A scratch script produced this
table (`vent_check.py`, not committed). It loads the cached BEM fields
through the project's own functions and does not change any score:

| Frequency | Cone alone: front − back | Cone + vents: front − back |
|---|---|---|
| 150 Hz | 3.3 dB | 6.0 dB |
| 200 Hz | 4.2 dB | 5.8 dB |
| 250 Hz | 4.7 dB | 5.6 dB |
| 300 Hz | 5.0 dB | 5.3 dB |
| 400 Hz | 7.1 dB | 5.5 dB |
| 500 Hz | 9.5 dB | 7.0 dB |

This suggests that most of today's rear rejection comes from the cabinet
shading the cone. The vents add 2.7 dB at 150 Hz, almost nothing at 300 Hz,
and reduce it from about 350 Hz up. Above about 180 Hz the vent term raises
the level behind the speaker instead of cancelling it. That suggests the
vents' timing at the rear is wrong. This is more likely a position or path
problem than a felt problem. It is a question for plan 2, not a settled
finding. (This split is not in the report.)

### 2. The bass needs about +14 dB of gain to match the coax

The crossover matches the bass to the coax LF with +14.2 dB of gain. The
arithmetic checks out. The datasheet sensitivities, plus the bass radiating
into full space while the coax is closer to half space, explain most of it.
But part of the 14 dB may come from the suspect 550–650 Hz coax peak (open
question 1), or from the de-embedding double-count.

- Each Dayton RSS265HO-44 is 84.5 dB at 2.83 V and 1 m (datasheet). The pair
  in parallel is about 6 dB more.
- The BMS coax is 96 dB.
- The bass works at 80–150 Hz, where the narrow cabinet radiates into full
  space. The coax is matched at 300–1000 Hz, where it already gets most of
  its baffle gain. That is several more dB.

In a passive design the coax would be padded down by about 14 dB. With
active drive, the bass amplifier channel needs that much more gain, or the
coax channels run 14 dB down. If the coax peak is overstated, the true gap
is a few dB smaller.

### 3. The puck model is already inside the trusted band

Above the trusted limit the coax is modelled as a "puck": the coax on a
round cylinder, with no cabinet facets or edges. Between half the trusted
limit and the limit, the cabinet result fades into the puck result. That
fade runs from **1.43 kHz to 2.85 kHz**, which is inside the trusted band.

So the trusted scores near the top of the band are partly model, not fully
checked 3D:

- presence_lift (trusted, 2–2.85 kHz): every point uses the puck; mean
  weight 0.72.
- di_smooth_lfhf (700–2800 Hz): 48% of points; mean weight 0.24.
- flat23, window_spread and rolloff_monotonic (150–2851 Hz): 24% of points;
  mean weight 0.13.

**Your choice before plan 2:**

- (a) Move the trusted limit down to where the fade starts (about
  1.4–1.9 kHz), so "trusted" means pure checked 3D.
- (b) Narrow the fade, for example to the top third of an octave below the
  limit, so less of the trusted band is puck.
- (c) Keep it as it is, and treat the disclosed puck share as part of the
  trusted result.

This matters most for the 2–5 kHz presence check, which you care about. Its
trusted part is already 72% puck.

## Other things the plots show

**Checks from the plan.**

- **The H polar map narrows with frequency.** Yes. Below about 180 Hz the
  whole speaker is within 6 dB of on-axis in every direction. The −6 dB
  angle is then about ±95° at 190 Hz, ±80° at 300 Hz, ±65° at 500–750 Hz,
  ±60° near 900 Hz and ±45° at 1.35 kHz (read from the cached fields, both
  sides the same).
- **The rear rejection plot shows the cardioid band.** Yes. The 150–300 Hz
  and 300–500 Hz bands are shaded. Read off `rear_rejection.png` (the
  150–500 Hz values agree with the cached fields in the table under
  observation 1): the 3D curve sits at about 5–6 dB from
  50 Hz to 450 Hz, then rises to about 18 dB at 850 Hz as the cone beams.
  The lumped curve is 20–25 dB below 300 Hz and falls to −1.4 dB at 900 Hz.
- **The spinorama sound power sits below the listening window.** Yes, at
  every frequency. Below 100 Hz they meet, as they should. The DIs rise
  smoothly from 0 dB at 50 Hz to about 5–9 dB at 2 kHz.

**Open questions.** These look odd. I did not change anything to chase them.

1. **A large coax LF peak at 550–650 Hz.** On-axis, listening window, early
   reflections and sound power all peak there (about 98 dB at 2 m) and fall
   about 8 dB by 1 kHz. Sound power follows, so it is source strength, not
   beaming. From the cached fields, the coax cone's on-axis pressure per unit
   volume velocity is +6 dB above a free-space point source at 300 Hz, and
   +13.5 dB at 570 Hz. A cone in a flat baffle can give about +6 dB, plus a
   dB or two of edge ripple. The datasheet (flat baffle) rises only about
   1 dB from 300 to 600 Hz. Two possible causes:
   - real: the recessed cone and its 290 mm flared mouth act as a short horn
     with a front-cavity resonance; or
   - model: the report already warns that de-embedding the datasheet with a
     flat piston may double-count the cone and cavity shape that the 3D
     model adds again.

   This also feeds the coax LF matching band (300–1000 Hz) and so the bass
   gain in observation 2.
2. **A deep on-axis notch in the coax LF field at 1.8–1.9 kHz.** The cone's
   3D field on axis drops about 15 dB there. The datasheet shows only a
   1–2 dB dip near 1.7 kHz. It barely shows in the total, because the LF way
   is already 12 dB down at 1.85 kHz after the 1.4 kHz crossover. It rests on
   the guessed cone depth (65 mm) and the cavity shape.
3. **At 220–450 Hz the deepest rear null is near ±120°, not 180°.** The H
   polar map shows a −10 dB island around ±120°. At 300 Hz the whole speaker
   is −12.3 dB at 120° but only −4.9 dB at 180° (from the cached fields).
   Straight behind, the sound coming round both sides of the cabinet
   probably arrives in step and fills the null. The gate is scored at 180°,
   which is the worst angle here. Assuming from your room description that
   the wall is straight behind the speaker, 180° is the right place to
   score. This just explains why rear_hemisphere and cardioid_hard both
   land at 180°.
4. **Below 100 Hz the speaker is about 2 dB louder below the axis than on
   it** (V polar, −30° to −120°). The bass drivers sit at 200 mm and 520 mm,
   well below the coax at 950 mm, so at 2 m the downward points are closer to
   them. This is plausible geometry, not an error.
5. **Above 2.85 kHz (low confidence, puck only):**
   - Around 3.5–4 kHz the ±30° directions are up to 5 dB louder than on
     axis. The window family shows a −5.5 dB on-axis dip at about 3.6 kHz.
   - The 23° curve has a single deep notch at about 11.5 kHz (about 15 dB).
   - A round puck makes every edge path arrive in step on the axis, so it
     tends to overstate on-axis peaks and notches. Treat these as shapes to
     re-check after the cone is measured, not as facts.
   - Above 12 kHz each direction keeps its 12 kHz pattern (the "puck hold"),
     so the polar maps show flat horizontal stripes and the DIs go flat. The
     two rolloff_monotonic warnings (12.7 and 16 kHz) sit in this held
     range.
6. **Is 2 m the right distance?** The scores use 2 m (spec §6.1): close to
   the 3 m seat, and common practice. But the wall that the cardioid is
   meant to protect is only 75 mm behind the back panel. At 75 mm the vents
   and the cone are in each other's near field, and the rear sound that
   reaches the wall is not the 180° point at 2 m. Should the rear rejection
   also be scored closer to the back panel, or on the wall plane?

## Provisional assumptions in force

- **Coax cone depth: 65 mm.** An estimate until the driver is measured with a
  profile gauge.
- **Coax HF exit: 25.4 mm.** Also an estimate. Above the trusted limit the
  coax's horn behaviour rests on these two numbers.
- **Coax cutout: 230 mm.** An estimate.
- **Fixed crossover:** acoustic LR4 at 170 Hz (bass to coax) and 1400 Hz
  (coax LF to HF). No EQ, no delays. Gains as in the table above.
- **Datasheet de-embedding:** the coax LF volume velocity above 300–500 Hz is
  the datasheet on-axis response divided by that of a flat piston (202 mm)
  in an infinite baffle at 1 m. The datasheet already contains the real cone
  and cavity, which the 3D model adds again, so this may double-count them
  (see open question 1). Below 300 Hz the lumped model is used. The two
  agree within 0.5 dB at 400 Hz.
- **Filler and coax volumes:** from the coax chamber design (spec §7.1.1).
  The coax is a basket cylinder until it is measured. The chamber volume,
  the lumped cardioid model and the vent strengths rest on these.
- **Vent resistance:** R = 190 Pa·s/m, shown with a ±30% spread.
- **Distance:** all scores at 2 m (see open question 6). The 50 m cardioid
  row is for reference.

## Run times and element counts

**Cabinet bands.** The full run reused these from the mesh check's cache, so
its log shows 0.0 s for them. The times below are from the mesh check log
(`out/mesh_check.log`) and `mesh_check.json`. "Full-run mesh" is the mesh
the scores use.

| Band ≤ | Freqs | Full-run mesh: h / nodes / elements / solve | Check mesh: h / nodes / elements / solve | Band total |
|---|---|---|---|---|
| 250 Hz | 29 | 50.0 mm / 2668 / 5332 / 367 s | 37.5 mm / 4423 / 8842 / 868 s | 1236 s |
| 500 Hz | 14 | 50.0 mm / 2668 / 5332 / 174 s | 37.5 mm / 4423 / 8842 / 418 s | 592 s |
| 1000 Hz | 12 | 50.0 mm / 2668 / 5332 / 149 s | 37.5 mm / 4423 / 8842 / 358 s | 508 s |
| 2000 Hz | 13 | 28.6 mm / 4132 / 8260 / 347 s | 21.4 mm / 7104 / 14204 / 962 s | 1310 s |
| 3000 Hz | 7 | 19.1 mm / 6966 / 13928 / 499 s | 14.3 mm / 12075 / 24146 / 1494 s | 1994 s |
| 4000 Hz | 5 | 14.3 mm / 10973 / 21942 / 893 s (full run only) | not run: needs about 24.0 GB, budget 22.0 GB | 893 s |

The 4000 Hz row's elements and time are from `out/full_run.log`. Its
frequency count, h and node count are from the cached field's `meta`, as
shown in the Run table of `out/directivity/report.html`.

Nothing in the report uses the 4000 Hz cabinet band. All its frequencies are
above the trusted limit, where the coax comes from the puck and the bass and
vents are not used. It costs about 15 minutes per full run.

**Puck bands.** Elements and solve times are from `out/full_run.log`.
Frequency counts, h and node counts are from the cached fields' `meta`, as
shown in the Run table of `out/directivity/report.html`.

| Band ≤ | Freqs | h | Nodes | Elements | Solve |
|---|---|---|---|---|---|
| 500 Hz | 1 | 30.0 mm | 1823 | 3642 | 7 s |
| 1000 Hz | 12 | 30.0 mm | 1823 | 3642 | 73 s |
| 2000 Hz | 13 | 28.6 mm | 1845 | 3686 | 85 s |
| 3000 Hz | 7 | 19.1 mm | 2171 | 4338 | 61 s |
| 4000 Hz | 5 | 14.3 mm | 2666 | 5328 | 66 s |
| 8000 Hz | 13 | 7.1 mm | 7297 | 14590 | 1001 s |
| 12000 Hz | 8 | 4.8 mm | 16345 | 32686 | 3168 s |

**Wall-clock** (from `out/round0_driver.log`):

| Step | Start | End | Time |
|---|---|---|---|
| Mesh check | 17:22 | 18:56 | 1 h 34 min |
| Full run | 18:56 | 20:25 | 1 h 29 min (report: 5357 s) |
| Total | | | 3 h 03 min |

The BEM test suite (`uv run pytest -m bem`, 12 tests) passed before the runs
in 146 s.
