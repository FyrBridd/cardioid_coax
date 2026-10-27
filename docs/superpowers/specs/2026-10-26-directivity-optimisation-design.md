# Directivity Optimisation Loop — Design

Date: 2026-10-26
Status: Draft for approval
Parent spec: `2026-10-22-cardioid-coax-floorstander-design.md` (the "main
spec"). This document replaces the main spec's planned "Plan 2" and covers its
§7.2 crossover targets and §7.3 items 2, 3, 4 and 6.

## 1. Intent

Get the best directivity this speaker can have, by changing its shape and
settings in rounds:

- the baffle and cabinet shape;
- the coax position and its flare;
- the cardioid vents;
- the crossover points.

**Why.** The owner listens 23° off-axis with no toe-in, in a room with nearby
walls. What they hear is mostly decided by directivity, not by the on-axis
curve. The shape of the box has to be frozen before the CNC cut, so this is
the moment to get it right.

**How the owner wants to work** *(the owner's words, summarised)*:

- The work is iterative. Whether a round runs in VituixCAD or in our own code
  doesn't matter.
- **The builder presents its best design first.** The owner then reviews the
  directivity and suggests changes. The owner's suggestions become the
  variants for the next round.

## 2. What one round is

1. **Variants.** The builder generates a set of design variants by changing
   the parameters in §5.
2. **Simulate.** Each variant is simulated: off-axis responses of every source
   (coax LF, coax HF, cardioid vents, each bass driver) on the real cabinet.
3. **Crossover.** Each variant gets a provisional crossover (§6.3), so it is
   judged as a complete speaker, not as bare drivers.
4. **Score.** Each variant is scored against the targets in §4.
5. **Report.** The builder writes a review report (§7) with **one
   recommended design**. The report says why that design beat the runners-up
   and what each change bought, in dB.
6. **Owner review.** The owner reads the plots and asks questions or suggests
   changes. Those suggestions are the starting variants of the next round.

Round 0 is the current design (`params/default.toml` at commit `94cffea`),
simulated and scored as the baseline. Every later report compares against
round 0 and against the previous round's best.

When the owner accepts a recommended design, it is written into
`params/default.toml`. The CAD, DXF and STEP files then regenerate from it as
they do now. The main spec is updated with the decided values.

## 3. Goals, non-goals, anti-goals

### Goals

- **A recommended geometry** that scores best on the targets in §4, within
  the constraints in §5, together with the evidence for it.
- **A loop that is cheap to repeat.** An owner's suggestion should become a
  scored variant without hand-work. A round should run unattended *(my call:
  overnight, at most about 12 hours on the owner's machine)*.
- **Honest confidence.** Every plotted curve and score shows where it relies
  on unmeasured driver behaviour (§8).
- **VituixCAD files** for the recommended design, so the owner can open it,
  check it and tune the passive crossover.
- **Active-mode crossover targets** as data: points, slopes, delays and EQ
  (main spec §7.2).

### Non-goals

- Passive crossover component values. The owner designs those in VituixCAD,
  from the exports (main spec §8).
- Panel resonances and bracing, fastening, the driver access route and the
  terminal panel. These stay with the later cabinet-construction plan.
- Room correction. The room model is used only to show room effects (for
  example the seat-behind-wall dip) so that they aren't blamed on the speaker.
- A final tuning. That needs the drivers and a measurement mic (main spec §9).

### Anti-goals

These would make it a failure even if the scores looked good:

- **A score that hides what the owner would hear.** A single number that
  rewards one thing while the 23° response gains 2–5 kHz lift, or that
  averages away a narrow problem. The report always shows the curves behind
  every score.
- **Confidence the simulation hasn't earned.** A geometry chosen because of
  detail above the 3D limit (2.85 kHz), where the coax's horn behaviour is
  guessed (§8), presented as if it were solid.
- **Breaking the main spec's anti-goals** in search of a better score:
  - visible driver frames;
  - a grille;
  - rounded top or bottom edges;
  - a cabinet that doesn't read as a floorstander;
  - a design that only works toed-in.
- **Losing the cardioid.** No variant may trade away the hard cardioid
  requirement for better directivity higher up.

## 4. Targets and scoring

The targets are the main spec's. The builder turns them into scores *(the
exact formulas are my call; the builder may refine them, and the report shows
the formula next to each score)*:

| Target (main spec) | Score |
|---|---|
| Flat at 23° horizontal | Deviation of the 23° response from a flat line, after the provisional crossover, over 150 Hz–16 kHz. Reported as max and RMS, in 1/3-octave bands. |
| Within ±1 dB of 23° over 0–30° horizontal | For each horizontal angle from 0° to 30° in 5° steps, the worst deviation from the 23° curve, in 1/6-octave smoothing. Reported per band, with the band where ±1 dB holds. |
| No 2–5 kHz elevation | For 0–30° horizontal and ±10° vertical: the largest amount any curve rises **above** the 23° curve's 2–5 kHz average, or above its own 1–2 kHz level. Any lift is flagged. Flat to slightly relaxed is acceptable. |
| Smooth directivity through both crossovers | The sound-power DI and the early-reflections DI, with the listening window re-centred on 23°. Score: the largest deviation of each DI from a smooth fit, within one octave of each crossover point. |
| Smooth roll-off beyond 30° | The horizontal response at 40–90° is monotonic with angle in each 1/3-octave band. Bands that break this are flagged. |
| Cardioid, hard | At least 10 dB rear (180°) rejection, 150–300 Hz. Pass/fail. |
| Cardioid, aim | At least 6 dB, 300–500 Hz. Reported, with a warning if missed. |
| Rear hemisphere | The worst rejection over 120–180° horizontal, 150–500 Hz (main spec: side walls are near). Reported. |

**Ranking** *(my call)*:

1. Variants that fail the hard cardioid requirement or a constraint in §5 are
   out.
2. The rest are ranked on the ±1 dB window and the 2–5 kHz check first,
   because these are the owner's stated priorities.
3. DI smoothness and 23° flatness break ties.

The ranking only decides the order. The recommendation is the builder's
judgement, argued in the report from the curves.

**Cardioid authority** *(my call)*: the main spec measures the hard
requirement with the lumped model. Once the 3D simulation includes the vents
(§6.1), its result is the one that decides pass/fail. The lumped model is
still used to tune the vent resistance. Where the two disagree by more than
3 dB, the report shows both.

**What the cardioid gate measures** *(from main spec §7.1 and goal 4)*:

- The **coax section alone**: the coax cone plus its vents, without the bass
  drivers and before the bass-to-coax crossover filter. The bass drivers
  radiate equally front and back, so the whole speaker can't reach 10 dB near
  the crossover. That is expected, and it isn't what the requirement means.
- The bands stay fixed at 150–300 Hz (hard) and 300–500 Hz (aim), wherever
  the crossover point moves within its range.
- *(My call.)* The whole-speaker rear rejection after the crossover is also
  plotted, for information only. It shows how the crossover choice affects
  the back-wall reflection.

**Scoring band and confidence** *(my call)*:

- **Ranking** uses only the band where the 3D results are trusted (§8),
  2.85 kHz and below (measured 2026-10-26, §8).
- **Above that limit**, the scores (including most of the 2–5 kHz check) are
  computed and shown in the report with the low-confidence marker. They don't
  rank variants.
- **Why the solver isn't pushed to 5 kHz.** Above the 3D limit the coax's
  own cone and horn decide the result. Their profile is guessed until the
  driver is measured, so a finer 3D mesh would add cost without adding
  trust.
- **Re-check after measurement.** The 2–5 kHz check is re-run when the cone
  profile has been measured with a profile gauge (main spec §9), before the
  CNC release. If a flare or baffle decision depends on that band, the
  report says so, and that decision stays provisional until the re-check.

## 5. What can change, and the constraints

### Constraints (fixed)

- **External width ≤ 330 mm**, **external depth ≤ 350 mm** *(owner,
  2026-10-26)*. The look may change within these.
- **Coax acoustic centre at seated ear height** (0.95 m for now; a parameter).
  Overall height follows from it.
- All anti-goals of the main spec (§3 above).
- **Rear-mounted drivers**, with the coax flare as a printed insert. The
  print service limits from the fillers work still apply: ABS up to
  580 × 480 × 480 mm, other materials up to 250 mm cubes.
- **Bass volume** at least 25 L net per driver.
- **The cardioid chamber** keeps the decided §7.1.1 scheme: printed fillers,
  side slots and a tuned resistance. Its volume, slot count and slot position
  may be re-tuned within that scheme. The hard cardioid requirement must
  still pass.
- **2.5D CNC panels.** Any curved or tapered shape must be buildable as flat
  routed panels plus printed or laminated parts. For example: a roundover cut
  with a router, a chamfer, or a printed edge trim or facet.

### Parameters the search may change

- **External width**, from about 305 mm to 330 mm. The bass frames are 266 mm
  across and mount from behind inside 18 mm side panels, which sets the lower
  limit. The builder computes the exact limit from the driver data and panel
  thickness.
- **Shape families the owner accepts** *(owner, 2026-10-26)*. The front face
  can only be narrower than the box through these:
  - **Vertical side-edge roundovers**, of any size the build allows.
  - **Large side chamfers:** up to about 30–45 mm at 45°, the full height of
    the vertical front edges. They are built from solid or printed edge
    strips. The box stays a plain rectangle seen from the side.
  - **Angled upper facets:** the baffle on each side of the coax angles back,
    like a shallow wedge. The bass section below stays flat-fronted. The
    facet angle, the facet width and how it meets the flat bass section are
    search parameters. The builder chooses how it is built, for example
    printed or laminated facet pieces on a flat routed baffle.
  - A plain flat front, with none of the above, is always one of the
    variants.
- **Rejected by the owner:** a coax section narrower than the bass section
  (a stepped two-part look).
- Top and bottom edges stay square in every family.
- **Coax horizontal offset** from the centre line, mirrored between the two
  cabinets. *(My call to include it: an offset spreads the edge-diffraction
  arrival times.)*
- **Coax distance to the cabinet top.**
- **Flare:** mouth diameter, depth (including a thicker outer baffle layer;
  this resolves the shallow-flare open item) and profile.
- **Cardioid vents:** position from the front, slot count and resistance,
  with the chamber volume re-solved by the fillers.
- **Bass driver heights** and spacing from the coax, within the bass volume
  constraint.
- **Crossover points:** bass-to-coax 150–250 Hz, and coax LF-to-HF
  1.0–2.0 kHz.

Every parameter listed lives in `params/*.toml`, so a variant is just a
params file. Parameters that don't exist yet are added there.

## 6. How it works

### 6.1 Simulation

- **3D boundary-element simulation** of the whole cabinet exterior, in free
  space *(my call: free space for scoring; the floor and room are applied
  afterwards only for the in-room plots)*.
  - Each source is simulated on its own, at unit strength:
    - the coax cone with the flare;
    - the HF exit;
    - the vent slots on each side;
    - each bass driver.
  - The lumped models (cardioid chamber, sealed bass, the drivers' electrical
    side) then set each source's real strength. The sound field is the sum.
  - So changing the vent resistance, the drive voltage or the crossover
    needs no new 3D run. Only a change of shape does.
- **Frequency range and angles:**
  - 3D run from 50 Hz to the highest frequency the mesh supports on this
    machine. Measured at 2.85 kHz on 2026-10-26 (16 cores, 30 GB RAM, no
    GPU), capped by memory: the 4 kHz check mesh needs about 24 GB (§8).
    Mesh density and symmetry tricks are the builder's choice.
  - Horizontal and vertical planes, 0–360° in 5° steps. Plus enough of the
    sphere for sound power (main spec §7.3 asks for 10° or finer).
  - At 2 m, and extrapolated to far field. *(My call: the scores use 2 m,
    which is near the 3 m seat distance and matches common measurement
    practice.)*
- **Above the 3D limit:** the coax's own directivity takes over. The builder
  models the cone and flare as a horn from the estimated cone profile. For
  example: an axisymmetric boundary-element model to 20 kHz, joined to the
  3D results with a crossfade over an octave. The method is the builder's
  choice. These results carry the low-confidence marker (§8).
  - **The crossfade starts at the 3D limit (owner's decision, 2026-10-27).**
    Trusted results come from the 3D model alone. The crossfade runs from the
    limit to an octave above it, or to the top of the 3D solve if that comes
    first, and it is at least a third of an octave wide. With the 3D solve
    ending near 3.8 kHz, it covers 2.85–3.8 kHz. The 3D results in the
    crossfade are not mesh-checked. The axisymmetric model has no vents, so
    the vents fade out across the crossfade.
- **Driver data:**
  - The datasheet on-axis and impedance curves of both drivers are digitised
    once into `data/` with their source noted.
  - The simulation supplies the cabinet's effect. It is combined with the
    datasheet response after allowing for the datasheet's own measuring
    conditions, as best they can be inferred. The builder states the
    assumption.
- **Validation before use** *(required)*:
  - The 3D solver is checked against cases with known answers: a piston in
    a sphere (analytic), and the baffle step of a simple box.
  - Mesh convergence is checked on the real cabinet: the scores must not
    move by more than 0.5 dB when the mesh is refined.
  - The report states the upper frequency where the 3D results are trusted.

### 6.2 Search

- **The search method is the builder's choice.** It has to fit the round
  budget (§3): 3D runs are the expensive step.
- *(My suggestion)*:
  - First screen many variants coarsely, at low mesh density and up to about
    1.5 kHz, where the baffle and cardioid matter most.
  - Then run the most promising few at full resolution.
  - Change several parameters at once, rather than one at a time, so that
    interactions show.
- **Each variant is reproducible:** its params file, its mesh settings and
  the code commit are stored with its results. Results are cached by those
  three, so re-running a round doesn't redo identical simulations.

### 6.3 Provisional crossover

So that geometry is judged fairly, each variant gets a crossover fitted to
its own simulated responses *(my call)*:

- **The idealised fully-active mode** is used for scoring: each of the three
  ways gets its own filters.
  - Acoustic target slopes: Linkwitz-Riley 4th order, as a starting point.
  - Plus a delay per way.
  - Plus parametric EQ fitted to make the 23° response flat.
  - The crossover points are searched within the ranges in §5.
- This is what the active-mode deliverable will be for the chosen design.
  It is exported as data: points, slopes, delays, EQ biquads in a plain
  table.
- **Passive and semi-active modes** are not scored in the loop. The report
  does flag where a passive network would struggle to copy the active one,
  for example EQ needing more than about 6 dB of boost, or very steep slopes.
  The owner then designs the passive network in VituixCAD from the exports.

### 6.4 VituixCAD exports

For the recommended design, and on request for any variant:

- **FRD files** for each way, for horizontal and vertical 0–180° in 5° steps.
  The file names follow VituixCAD's off-axis naming, so that VituixCAD
  imports the whole set at once.
  - The **coax LF set includes the vent radiation**, because the vents are
    driven by the cone's rear. Otherwise VituixCAD would show the coax
    without its cardioid.
  - Bass: each driver separately, and the pair.
- **ZMA files** for each way, from the lumped models, including the effect
  of the cardioid chamber on the coax LF impedance.
- Each file carries a header line marking it as simulated, with the
  confidence band.
- *(My call)* The FRDs are at 2 m in free space. A short note says how to
  set VituixCAD's distance and baffle options so it doesn't apply its own
  baffle simulation on top.

## 7. The review report

Each round produces one report *(my call: a local HTML page with PNG plots,
in `out/rounds/<n>/`; it can also be published as a private web page to view
on another device, if the owner asks)*. It contains, in this order:

1. **The recommendation.** The recommended design's key dimensions, and a
   picture of it next to round 0. A short argument for why it beat the
   runners-up, with the dB each change bought.
2. **The curves for the recommended design:**
   - horizontal and vertical polar maps;
   - the spinorama with the listening window re-centred on 23°;
   - the 0–30° horizontal family normalised to 23°, with the ±1 dB band and
     the 2–5 kHz band marked;
   - rear rejection vs frequency, with the hard and aim lines;
   - both DIs.
3. **Score table** for the recommendation, the runners-up, round 0 and the
   previous round's best. Pass/fail and warnings are shown.
4. **Confidence:** where each curve moves from trusted to estimated.
5. **Crossover:** the provisional points, slopes, delays and EQ, and any
   passive-mode warnings.
6. **What I'd try next**, so the owner can agree, redirect or add
   suggestions.

The owner's feedback can be given in chat. The builder turns it into the next
round's variants.

**What runs unattended** *(my call)*:

- The overnight job produces the simulations, scores, plots and a draft
  report. The draft has the score table and a ranking.
- The written recommendation and "what I'd try next" are written in a
  following session, from those results.
- Turning the owner's suggestions into variants is also done in a session.
  No rule-based write-up or suggestion parser is needed.

## 8. Confidence

Marked on every curve and score:

- **Trusted:** below the 3D limit (2.85 kHz, measured 2026-10-26 by the mesh
  check; capped by memory, since the 4 kHz check mesh needs about 24 GB
  against a 22 GB budget), for the shape effects: baffle, edges, cardioid,
  crossover-region lobing.
- **Estimated:**
  - The cardioid vent resistance as built. The felt is unknown, so the
    report shows a ±30% resistance band *(my call)*.
  - The coax displacement and filler volumes (a basket cylinder, until the
    driver is measured).
- **Low confidence:**
  - Above the 3D limit. The coax's horn behaviour comes from an estimated
    cone profile. BMS publishes no coverage data.
  - Datasheet curves were measured under conditions that are only partly
    known.

A geometry decision must not rest on low-confidence detail alone. If two
variants differ only there, the report says so and prefers the simpler build.

## 9. Work split

Three plans, built in order. Each produces something usable on its own:

1. **Simulation and scoring.**
   - Driver data digitised.
   - 3D solver and meshing from the existing CAD geometry, validated (§6.1).
   - Coupling to the lumped models.
   - The HF model.
   - Spinorama and the §4 scores.
   - Round 0 scored and reported in a first, simple report.
2. **Search and review.**
   - New parameters (§5) in params and CAD.
   - Variant generation, caching, the search.
   - The full review report (§7).
   - The first real round and recommendation.
3. **Crossover and exports.**
   - Provisional active crossover fitting (§6.3) inside the loop.
   - Active-mode data export.
   - VituixCAD FRD/ZMA exports (§6.4).
   - *(My call.)* **An early import check.** This machine has no Wine, so
     VituixCAD can't be run here. At the start of plan 3, a small test set
     goes to the owner to import on Windows: a few angles of one way, plus
     one ZMA. That confirms the naming and headers before the full exporter
     is built.

*(My call on order.)* Plan 2 needs a crossover to score fairly. Until plan 3
lands, plan 2 uses a fixed ideal crossover at the main spec's starting points
(170 Hz and 1.4 kHz, LR4, no EQ). Plan 3 then replaces it, and the round is
re-scored.

## 10. Left to the builder

- The 3D solver, mesher and HF model libraries. They must be open source and
  run on Linux. *(My suggestion: bempp-cl with its CPU backend, and gmsh for
  meshing.)*
- Mesh density, symmetry use, frequency spacing and the 3D upper limit, all
  justified by the validation in §6.1.
- The search method, within the round budget.
- Code layout, following the existing `src/speaker/` structure.
- Refinements to the score formulas in §4, provided the report shows the
  formula and the curves.
- How the datasheet's measuring conditions are allowed for (§6.1), stated in
  the report.
