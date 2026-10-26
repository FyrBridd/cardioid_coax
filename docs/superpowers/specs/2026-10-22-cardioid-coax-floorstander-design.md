# Cardioid Coax Floorstander — Design

Date: 2026-10-22
Status: Draft for approval

## 1. Intent

A DIY floorstanding loudspeaker in the spirit of the Fyne Audio F703SP: a 10"
compression-driver coaxial for mids and highs, with two sealed 10" bass drivers
below it. The design is done with AI assistance as a parametric, code-driven
pipeline. The output is a 3D model and simulated measurements, good enough to
send the cabinet for CNC cutting and to design crossovers in VituixCAD.

**Why this design:** the speakers live in a large open-plan living room with a
hard floor and ceiling. They sit close to the back wall and near some side
walls, and the room tends to "suck the bass out". The owner wants:

- **Constant directivity**, so the sound is even over a reasonably wide
  listening area and the room reflections sound like the direct sound.
- **Passive cardioid loading of the coax** through the lower midrange, to cut
  the back-wall reflection (SBIR) and room excitation from speakers placed near
  the wall.
- **Sealed bass**, placed and tuned to excite room modes as little as possible
  while recovering the low end the room loses.

## 2. Goals

1. **Neutral tonality at the real listening position.** Speakers are 2.5 m
   apart, with the main seat about 3 m away. There is **no toe-in** (a WAF
   requirement), so the listener sits about **23° off-axis horizontally**.
   - The response is flat at 23° horizontal.
   - It stays within about ±1 dB of that from 0° to 30°, so nearby seats and
     the on-axis direction are never hotter in 2–5 kHz *(my call for the
     window)*.
   - There is a smooth, monotonic roll-off beyond 30°.
   - The spinorama listening window is redefined as 0–30° horizontal and ±10°
     vertical.
2. **No upper-mid or presence elevation.** The owner is very sensitive to
   2–5 kHz lift. In this band, flat to slightly relaxed is acceptable; elevated
   is not.
3. **Constant (smoothly varying) directivity** through the coax-to-bass and the
   coax LF-to-HF crossovers. No directivity flare at either crossover.
4. **Cardioid behaviour** of the coax section from the bass crossover (about
   150 Hz) up to where the baffle and cone become naturally directional (about
   600 Hz and above).
5. **Output:** clean at loud levels for metal, and composed for vocal music
   (London Grammar, Florence and the Machine). Eventual target in fully active
   mode: about 100–105 dB peaks at the seat for the pair. The first passive
   version is not expected to reach maximum output.
6. **Three drive modes**, selectable on the same cabinet, so the system can grow:
   - **Passive:** a single amp (the owner's JBL SA550) drives a full passive
     crossover.
   - **Semi-active:** the coax's internal LF/HF crossover is passive; the
     coax-to-bass crossover is active (DSP plus separate amp channels).
   - **Fully active:** coax HF, coax LF and bass are each amplified and crossed
     over in DSP.

## 3. Non-goals

- Designing amplifiers, DSP hardware or DSP firmware. The design supplies
  target filters and EQ for the active modes; implementing them on a specific
  DSP is out of scope.
- Room EQ or multi-sub integration.
- Cost optimisation. Drivers are chosen on merit. The owner will revisit cost
  later, and has already rejected the Purifi bass drivers as not worth it here.
- A final, measurement-verified tuning. That needs a measurement mic, which the
  owner doesn't have yet (see §9).

## 4. Anti-goals

These would make the result a failure even if it otherwise works:

- Audible upper-mid or presence elevation at the listening position.
- Visible driver mounting frames or plates on the front.
- Rounded top or bottom edges (rounded vertical side edges are fine).
- A cabinet that doesn't read as a floorstander, or as a very big standmount on
  a very short stand.
- A design that only sounds right on-axis (toed-in). It must work facing
  straight ahead.

## 5. Constraints

| Item | Constraint |
|---|---|
| External width | about 320 mm |
| External depth | **350 mm max** |
| Height | Set by the coax: **acoustic centre at seated ear height** on the owner's sofa. Working assumption is 0.95 m to the coax centre, giving about 1.1 m overall. Parametric, so it can be updated once the owner measures. |
| Front edges | Top and bottom edges square. Vertical side edges may be rounded or chamfered if it helps diffraction. |
| Driver frames | Hidden: rear-mounted, or front-mounted with a cover. Rear-mounting is preferred, and the coax cutout may double as a waveguide extension. |
| Fabrication | Flat panels cut by a CNC service and assembled by the owner with basic woodworking tools. 3D-printing services are available for parts such as the waveguide/flare insert, vent pieces and trims. |
| Amplification today | JBL SA550 (about 90 W into 8 Ω, about 150 W into 4 Ω). In passive mode the load must be safe for it, with a minimum of about 3.5 Ω. |
| Measurement | None yet. The design must be buildable from simulation alone, with a defined path to verify with a mic later. |
| Tools | Linux workstation. VituixCAD (Windows) for crossover design. The pipeline must otherwise use free/open-source tools. |

### 5.1 Room

These figures are estimated from the owner's floorplan (scale about
110 px/m on the plan), unless noted otherwise. All of them must be parameters
in the pipeline.

| Item | Value |
|---|---|
| Sitting room | about 4.6 m (along the speaker line, N–S) × 4.0 m (the listening axis, E–W). Open plan: the north wall has a wide opening into a 3.0 × 2.5 m dining room, and there are doors to the hall (SW corner) and the kitchen (south wall). |
| Ceiling | about 2.0 m, the owner's estimate ("pretty standard"). To be measured. |
| Floor / ceiling | Hard (reflective) |
| Speakers | Against the west wall (the stair wall), firing east. The plan shows the cabinet backs about 0.2–0.3 m from the wall. Spacing: the owner states 2.5 m. The plan sketch looks closer, so confirm on site. |
| Side boundaries | North speaker: a solid wall section about 0.8 m to its side; the dining-room opening starts beyond it. South speaker: a short wall return (the hall jog) about 0.4–0.5 m away, then the hall door. The two sides are **asymmetric**. |
| Seat | About 3 m from the speakers, centred between them, about **0.4–0.6 m from the east wall** behind the listener |

What follows from this *(my analysis, to be confirmed by the room/boundary
simulation in §7.3)*:

- **Front-wall SBIR:** the null falls at about c/(4·d), where d is the distance
  from the acoustic centre to the wall.
  - Cabinet back 0.25 m from the wall (d ≈ 0.6 m): **null at about 145 Hz**.
    That lands on the bass crossover, where the bass drivers are omni and the
    cardioid can't help. This is the worst case.
  - Cabinet back about 0.05–0.10 m from the wall (d ≈ 0.4–0.45 m): null at
    about 190–215 Hz. That is above the crossover, inside the cardioid band,
    where the rear rejection suppresses it.
  - **Recommended placement is therefore as close to the wall as practical,
    and the design should assume it** (cabinet back 50–100 mm from the wall).
    The rear terminal panel must allow this, so the terminals are recessed or
    angled.
  - The bass drivers then get boundary gain below the crossover, which helps
    the passive-mode level matching.
- **Seat-behind-wall reflection:** with the listener about 0.5 m from the east
  wall, there will be a dip around 150–200 Hz at the seat. The speaker can't fix
  it (room treatment or seat position only). The spinorama and in-room
  predictions should show it, so it isn't misattributed to the speaker.
- **Axial modes:**
  - E–W (4.0 m): 43 / 86 / 129 Hz. The seat is near the east wall, at a
    pressure maximum for all of them.
  - N–S (4.6 m, partly open): about 37 / 75 Hz.
  - Vertical (2.0 m): about 86 Hz.
- **Bass driver heights** are chosen to reduce coupling to the vertical mode
  (node at mid-height) and to spread the floor-bounce nulls of the two drivers.
- **Floor bounce at the seat** from the coax (ear height 0.95 m, 3 m away): null
  at about 310 Hz.

## 6. Decided by the owner

### 6.1 Drivers

**Coax: BMS 10C262 (8 Ω LF / 8 Ω HF)**
([datasheet](https://bmsspeakers.com/wp-content/uploads/2022/01/bms_10c262_datasheet.pdf)).
The owner chose it over the B&C 10CXN76 after reading the response plots: the
B&C's LF breakup runs about 17 dB from about 1.5 to 4 kHz, against the BMS's
milder 6 dB breakup from 2 to 4 kHz, and the BMS tweeter response is less
choppy.

| Parameter | Value |
|---|---|
| Fs / Qts / Qes / Qms | 65 Hz / 0.37 / 0.40 / 4.6 |
| Vas | 25.5 L |
| Xmax | ±5.5 mm |
| Sd | about 320 cm² |
| Re (LF) / Le | 5.9 Ω / 0.5 mH |
| Sensitivity | LF 96 dB, HF 112 dB (1 W/1 m). HF plot averages about 104 dB above 4 kHz. |
| HF | 1.75" polyester diaphragm. Recommended crossover 1.2 kHz. |
| Frame | OD 260 mm, about 10 mm flange, 146 mm overall depth, rear basket Ø228 mm, magnet Ø156 mm, 5.05 kg (ferrite). Cutout not stated; about 230 mm inferred from the drawing. |
| Coverage | **Not published.** No polar data. |

Known issues the design must handle:

- **HF on-axis hump** from 2 to 3.3 kHz, falling about 6 dB by 5 to 6 kHz.
- **LF cone breakup:** a dip at 1.6 kHz and a peak at about 2.7 kHz.
- **HF impedance resonance** at about 750 Hz, close to the 1.2 kHz crossover.
- **Limited LF excursion.** The coax must be crossed to the bass drivers at
  **about 150 Hz or higher**, because cardioid cancellation costs output at the
  bottom of its range.
- A forum builder reports a cone/compression-driver timing offset that needs
  delay compensation in active mode.

**Bass: 2 × Dayton Audio RSS265HO-44 per cabinet**
([product page](https://www.daytonaudio.com/product/1108/rss265ho-44-10-reference-ho-dvc-subwoofer-4-ohms-per-coil)),
front-firing and sealed.

| Parameter | Value |
|---|---|
| Fs / Qts / Vas | 27.6 Hz / 0.45 / 34.8 L |
| Xmax | ±13.25 mm |
| Sd | about 352 cm² |
| Voice coils | Dual 4 Ω. Coils are wired **in series** (8 Ω per driver, Re about 6.4 Ω, Le about 3.1 mH); the two drivers are wired **in parallel**, giving about 4 Ω for the pair. |
| Sensitivity | about 84.5 dB per driver (series coils), about 90.5 dB for the parallel pair (2.83 V) |
| Frame | OD 266 mm, cutout 235 mm, depth 136.5 mm |
| Sealed alignment | 25 L each: Fc 42.7 Hz, Qtc 0.70, F3 about 43 Hz. 30 L each: Fc 40.6 Hz, Qtc 0.66, F3 about 44 Hz. Before boundary/room gain and stuffing. |

Rejected: Purifi PTT10.0X08 (not worth about 4× the cost for use below
150 Hz), and Visaton TIW 250 XS (F3 about 72 Hz sealed in 25–30 L; it suits a
vented or larger box).

### 6.2 Layout and cosmetics

- **Coax at the top**, acoustic centre at seated ear height.
- **Bass drivers front-firing** on the same baffle below the coax, mounted from
  behind (hidden frames). Opposed side-firing was considered and rejected: the
  owner prefers a clean front with no side grilles.
- **Coax mounted from behind**, with the baffle cutout profiled as a flare that
  continues the cone's profile outward, acting as a shallow waveguide extension.
- Square top and bottom. Side edge treatment is set by diffraction simulation.
- **No grille.** The baffle, the flare and the visible driver cones are the
  finished face. No grille fixings are needed.
- **Cardioid vents may be visible on the side panels.** The owner prefers
  **several small slots** over one large opening (preferred, not critical).
  Both cabinets have the same vents on both sides, for a symmetric horizontal
  pattern.
- **The owner doesn't have the drivers yet.** Datasheet drawings are good
  enough to develop the design, but the final CNC release waits until the real
  drivers are measured (§9).

### 6.3 Deliverables

- A parametric 3D model of the full cabinet: baffle, coax flare, cardioid
  chamber and vents, bass chamber(s), bracing, damping placement, terminal
  panel.
- **STEP** files of the assembly and of each part, for 3D review and printing.
- **DXF** files of every flat panel for the CNC shop, with cutouts, rebates,
  dados and drill holes.
- Simulated measurements (see §7.3).
- Driver and system responses exported in VituixCAD-ready form (FRD/ZMA,
  off-axis sets), so the crossovers can be designed in VituixCAD.

## 7. The design in more detail

### 7.1 Cabinet architecture

- **Upper section: the coax and its cardioid chamber.** This is a separate
  enclosure behind the coax, sealed from the bass section, with **resistive
  vents**: openings loaded with acoustic resistance such as felt, foam or
  fabric over slots. Rear radiation passes through the vents with the right
  delay and level to cancel the front wave behind the speaker. The starting
  chamber volume is about 10 L (Qtc about 0.7 if it were sealed). The final
  volume, vent area, vent position and resistance come from simulation.
- **Vent position:** the cabinet back will be 50–100 mm from the wall (§5.1),
  so the vents go on the **side panels near the rear**, as an array of several
  small slots per side. Rear-panel vents would be loaded by the wall.
  *(Builder to confirm by simulation.)*
- **Cardioid target** *(revised 2026-10-26, agreed with the owner)*:
  - **Hard requirement:** at least 10 dB rear (180°) rejection from 150 to
    300 Hz, in the lumped model, for the chamber *as built*.
  - **Aim:** at least 6 dB from 300 to 500 Hz.
  - The original "10 dB to 600 Hz" target is not reachable with a side-vented
    10 L chamber in this cabinet. The aim above replaces it.
  - Rejection across the rear hemisphere (about ±120–180°) matters too,
    because of the nearby side walls, so a hypercardioid-leaning pattern is
    acceptable.
  - If the coax's ±5.5 mm Xmax can't support this at 150 Hz, move the
    bass-to-coax crossover up towards 200 Hz rather than weaken the rejection.
- **Coax chamber fillers and vents** are decided in §7.1.1.
- **Lower section: sealed bass.** Two Daytons, target 25–30 L net each
  (50–60 L total). Use as much of this volume as the cabinet allows *(my
  call)*, because it lowers Qtc and leaves room for the passive inductor's
  resistance.
- **Low-frequency target** *(my call)*:
  - Passive mode relies on the natural sealed roll-off (F3 about 43 Hz) plus
    the wall and floor gain.
  - The active modes may use a Linkwitz-transform style bass EQ. That is
    speaker EQ, not room EQ, so it is in scope. Target F3 is about 30–35 Hz,
    within excursion limits.
  - The 100–105 dB peak target at the seat applies down to 40 Hz in active
    mode. They can share one volume or have separate chambers. The
  two drivers' vertical positions are chosen to reduce floor-bounce cancellation
  and to stagger their coupling to floor and ceiling modes.
- **Baffle:** thick enough for rear-mounting plus a meaningful flare depth. The
  rear-mount cutout depth acts as a short tube in front of the cone, so its
  profile must be shaped (flared or chamfered), not a straight bore.
- **Bracing and damping** are set by panel-resonance simulation, with bracing
  that doesn't block the cardioid path or the bass drivers' rear volume.
- **Terminal panel** on the rear, supporting the three modes (§7.2).

### 7.1.1 Coax chamber: fillers and vents (decided 2026-10-26)

**Why.** Plan 1 showed that the upper section as drawn has about 18.8 L net
behind the coax. That gives only about 5 dB rear rejection over 150–300 Hz.
A lumped-model sweep compared the options. The owner chose to **keep the single
box and fill it down to 10 L with printed fillers**, together with doubling the
vent area. As built, the model gives 22.6 dB minimum over 150–300 Hz and
8.9 dB minimum over 300–500 Hz. With the doubled vents but no fillers it gives
13.5 dB and 4.3 dB, so the fillers mainly buy margin and the 300–500 Hz aim. Coax excursion stays at about 1.8 mm at 150 Hz and
25 Vrms. Rejected alternatives:
- a partition behind the coax, which clashes with the side vents;
- a shallower stepped upper section (best in the model, but it changes the look).

**Ranking of levers** (from the sweep, for anyone re-tuning): vent area has
the biggest effect. Next is placing the vents as far behind the cone as
possible. Smaller chamber volume comes third.

**Vents** *(decided)*:
- 12 slots per side, each 100 mm long (front-to-back) and 10 mm high, at a
  22 mm vertical pitch. Centred 270 mm from the front face and at the coax
  centre height, as now.
- Total open area: 24,000 mm² for both sides together.
- Slot fill resistance is re-tuned by the model. It is expected to be around
  190 Pa·s/m.
- *(My call.)* The tuner runs on every build:
  - It maximises the minimum rejection over **150–300 Hz only**, which is
    the hard requirement.
  - It never gives up 150–300 Hz margin to chase the 300–500 Hz aim. A miss
    on the aim is reported as a warning, not traded for.
  - The tuned value is reported in the summary, not written back to the
    params file.
- The default parameters change accordingly: `slots_per_side = 12`, with
  `specific_resistance_rayl` set to the tuned value.

**Fillers.** There are two kinds, both 3D-printed:

1. **Four corner pieces** *(decided)*:
   - One in each corner of the upper section's cross-section, outside a
     cylinder of the coax frame diameter plus clearance, coaxial with the coax.
   - The concave inner face keeps the full frame-diameter sweep clear, so the
     coax can be slid in from the back and bolted to the baffle with all
     fillers in place.
   - Each runs from the back face of the baffle to **15 mm in front of the
     vent slots**. With the current numbers, the slots start at y = 220 mm,
     so the pieces end at y = 205 mm.
   - The pieces don't extend into the vent band, so the air path from the
     back of the cone to every slot stays open.
   - They stay 40 mm clear of the coax's horizontal and vertical centre
     lines, so no part tapers into an unprintable sliver *(plan call)*.
   - Together they give 3.85 L.
2. **One rear block** *(decided)*:
   - Mounted on the inside of the back panel, behind the magnet. It comes out
     with the back panel.
   - It fills the rest of the volume needed to reach the 10 L net target
     (4.98 L, 234 × 127 × 167 mm, at the current numbers).
   - It is **sized by the model, not fixed**: its size is solved so that the
     chamber's net volume, measured from the geometry, equals
     `coax_chamber_net_volume_l`.
   - Clearances:
     - at least 10 mm behind the magnet;
     - at least 25 mm between the block and each side panel along the whole
       vent band, so the air path to the slots is not throttled.
   - If it can't be made big enough within those clearances, the build fails
     validation with a clear message. The volume is never quietly missed.

**Access** *(decided)*:
- The back panel over the coax chamber is **removable**, and is how the coax
  is installed and serviced.
- It must seal airtight when closed, for example with a gasket.
- Its fasteners and the gasket detail are left to Plan 3, together with the
  other fastening questions.

**Print and material** *(decided)*:
- The fillers are **ABS**, printed by the owner's print service. It offers
  build volumes up to 580 × 480 × 480 mm in ABS, or 250 mm cubes in other
  materials. ABS softens at about 95–100 °C, far above anything the chamber
  will see.
- The parts must be **airtight**, because air inside them would still count
  as chamber volume. Print with at least 4 perimeters, at least 5 solid
  layers top and bottom, and about 20% infill. Seal the surface by acetone
  vapour smoothing or by brushing on ABS slurry.
- At the current numbers every piece is under 250 mm in every direction, so
  no splitting is needed. If a re-tune makes a part larger than 250 mm, it
  can still be one ABS part up to the larger build volume.
- Fixing:
  - The corner pieces are bonded to the panels with PU adhesive.
  - The rear block is screwed to the back panel from outside.

**What the software must do** *(decided)*:
- Model the fillers as real CAD parts, and export them as STEP and STL for the
  print service. They carry the same `_PROVISIONAL` suffix as the other files
  until the drivers are measured, because their inner faces depend on the
  coax frame and magnet sizes.
- Compute the chamber's as-built net volume from the geometry. That is the
  upper-section interior less the coax displacement less the filler solids.
  Feed this volume to the cardioid model.
  - *(My call.)* Each filler subtracts its full outer envelope, because the
    parts are sealed.
  - *(My call.)* The rear block is solved to hit the target within ±0.05 L.
- The summary reports:
  - the as-built volume;
  - the tuned resistance;
  - the minimum rear rejection over 150–300 Hz and over 300–500 Hz.

  It warns if the result misses the hard requirement or the aim.
- This replaces Plan 1's ruling R9, under which the filler was not modelled.

### 7.2 Drive modes and crossover topology

| Mode | Amp channels per speaker | Coax HF/LF | Coax-to-bass |
|---|---|---|---|
| Passive | 1 | Passive network | Passive network (bass low-pass, coax high-pass) |
| Semi-active | 2 (coax, bass) | Passive network | DSP |
| Fully active | 3 (HF, LF, bass) | DSP | DSP |

- The passive coax network and the passive bass low-pass are **separate
  sub-networks**. In the active modes, link bars or jumpers on the terminals
  bypass the parts that aren't used. *(My call: links rather than a high-current
  switch, for reliability. The builder may propose otherwise.)*
- Passive-mode requirements:
  - HF padding of about 10–15 dB.
  - Contouring of the 2–3.3 kHz HF hump.
  - Handling of the 750 Hz HF resonance.
  - A low-pass on the coax LF steep enough to suppress the 2.7 kHz breakup.
  - Level matching of the coax LF (96 dB) to the bass pair (about 90.5 dB):
    lean on boundary/back-wall gain for the bass where possible, and pad the
    coax where not.
  - A minimum impedance safe for the SA550.
- Active-mode deliverables: target acoustic slopes, crossover points, delays
  and EQ curves, as data a DSP can be programmed from. No specific DSP is
  assumed.
- Crossover frequencies are to be optimised. Starting points: **bass-to-coax
  about 150–200 Hz**, and **coax LF-to-HF about 1.2–1.6 kHz**, the latter
  trading the 750 Hz HF resonance against directivity matching.

### 7.3 Simulation pipeline and outputs

The pipeline is code-driven and parametric: changing a dimension regenerates
the geometry and the simulations. *(My call: Python, using CadQuery or an
equivalent scripted CAD tool for geometry, and open-source solvers. The builder
picks the specific libraries.)*

Required simulated outputs:

1. **Bass alignment:**
   - sealed response, excursion and max SPL vs frequency for both drive levels
     (SA550 passive; active)
   - the effect of being close to the back wall and floor
2. **Coax on the baffle:**
   - horizontal and vertical off-axis responses from 0 to 180° in 10° steps or
     finer, including the rear-mount flare and baffle diffraction
   - effect of edge treatment
3. **Cardioid behaviour:**
   - front-to-back ratio and polar pattern vs frequency for the coax section
   - sensitivity to vent area, vent resistance and chamber volume
4. **System spinorama**, CTA-2034-style:
   - on-axis, listening window (re-centred on 23°), early reflections, sound
     power, and directivity indices
   - polar maps, horizontal and vertical
5. **Panel resonances:** modal frequencies and shapes of the main panels, with
   and without the proposed bracing.
6. **VituixCAD exports:** FRD/ZMA in VituixCAD's off-axis naming convention.
   - Coax responses come from the datasheet curves combined with the simulated
     baffle, flare and cardioid behaviour.
   - Bass responses come from the T/S parameters plus baffle effects.

The simulation must be honest about its limits. The cardioid vent resistance
and the damping are the least predictable parts, and the outputs should mark
which results are high-confidence and which need measurement.

## 8. Left to the builder

- Choice of specific open-source tools (CAD library, BEM/FEM solvers, lumped
  modelling) and the repo layout.
- Panel material and thickness. *(My suggestion: 18 mm birch ply or MDF, with a
  doubled or thicker front baffle for the rear-mount flare.)*
- The resistive material for the vents, chosen to match the tuned
  resistance. Chamber volume, vent geometry and vent location are now decided
  in §7.1.1.
- The exact cross-section shape of the filler pieces, within the clearances
  in §7.1.1.
- Bass chamber split (shared or separate) and the vertical positions of the
  two bass drivers.
- Bracing pattern and damping material placement.
- Edge treatment of the vertical side edges.
- Crossover topology details and component values. These are designed in
  VituixCAD by the owner, assisted by the pipeline's exports.
- CNC detail beyond the following *(my call, because no shop has been
  chosen)*:
  - Assume 2.5D flat-panel routing only: through-cuts, pockets, rebates and
    dados.
  - Inside corners have at least a 3 mm tool radius (dog-bones where needed).
  - DXF uses one layer per cut depth, with the depth in the layer name.
  - **The coax flare is a separate 3D-printed insert** fitted into a simple
    routed opening. That keeps the baffle 2.5D, allows correction once the
    real driver is measured, and gives a clean finished face with no grille.
- Terminal and link layout for the three modes.

## 9. Verification path

- **Before measurement:** simulation outputs only, with confidence marked.
- **When the drivers arrive (before CNC release):**
  - Measure the frame OD, flange thickness, bolt circle, basket and magnet
    clearances with calipers.
  - Take the coax cone and surround profile with a profile gauge.
  - Update the parametric model, and regenerate the cutouts, the flare insert
    and the DXFs.
  - Optionally sweep the impedance (for example with a DATS or a soundcard
    jig) to confirm the T/S parameters.
  - **DXFs are marked "provisional" until this step is done.**
- **Recommended purchase:** a calibrated USB measurement mic (e.g. miniDSP
  UMIK-1) before final tuning.
- **With a mic:** in REW, measure each driver in-cabinet, horizontal 0–90°
  (and rear) responses, and the cardioid front/back ratio. Feed the
  measurements back into VituixCAD, replacing the simulated FRDs, then refine
  the vent resistance and crossovers.

## 10. Open items to confirm with the owner

To be measured on site. None of these block starting work, because they are
parameters.

- Actual seated ear height on the sofa (currently assumed to be 0.95 m).
- Ceiling height (currently 2.0 m).
- Speaker spacing, as placed (stated as 2.5 m; the plan sketch suggests less).
- The achievable distance from the cabinet back to the wall (target
  50–100 mm).
