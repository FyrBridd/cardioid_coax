# Passive crossover: design

Date: 2026-10-30. Status: draft for the owner's approval.

Main spec: `docs/superpowers/specs/2026-10-22-cardioid-coax-floorstander-design.md` (§2 item 6, §7.2).
This amends the directivity spec (`2026-10-26-directivity-optimisation-design.md` §3), which left
the passive component values to the owner in VituixCAD. The owner asked on 2026-10-30 for as much as
possible to be done without VituixCAD, so the pipeline now designs the passive network itself.

## 1. Why

The first version the owner builds is passive, driven by one JBL SA550 channel per speaker (main spec
§2: "The first passive version is not expected to reach maximum output"). Only the fully active
crossover exists so far (round 3c: 150 Hz / 1.2 kHz). The passive network also decides how much room
the crossover boards need in the cabinet and how the terminal links are wired, so it comes before the
cabinet-construction work.

## 2. Goals

- **Component values** for the passive mode of the current design (`params/round-3.toml`), in
  standard values with realistic losses, and a parts list.
- **The same scores as the active crossover** (spec §4 of the directivity spec, scoring version 3),
  side by side with the active fit, so the owner sees what passive costs.
- **The semi-active mode** checked too: it uses the same passive coax network, with the bass-to-coax
  crossover in DSP (main spec §7.2). Its DSP settings come out as a sheet like `crossover.txt`.
- **Visual output the owner can read without VituixCAD:** schematics, the 23° response, polar maps,
  rear rejection, and the impedance against the SA550's limit.
- **A VituixCAD check later:** the network as a plain netlist the owner can type into VituixCAD.

## 3. Non-goals and anti-goals

Non-goals:
- A VituixCAD project file (.vxp). Its format is not documented well enough to write blind.
- Choosing a shop or a brand of parts. The parts list gives value, type, tolerance, DCR and power
  rating.
- Final tuning. That still needs the drivers and a mic (main spec §9).

Anti-goals (a failure even if the scores look good):
- **An unsafe load for the SA550.** The minimum impedance of the whole speaker must stay at or above
  3.5 Ω (main spec §7.2) from 20 Hz to 20 kHz, including the drop where the bass and coax branches
  overlap near 150 Hz.
- **A network that only works with exact values.** It must hold up with ±5 % capacitors and inductors.
- **Fixing problems above the trusted limit (2851 Hz) with narrow filters.** The data there is
  low-confidence, as the above-limit active PEQs showed.
- **Hiding the cardioid cost.** A passive network cannot delay the bass, so the whole-system rear
  rejection (9.9 dB active) will change. It is scored and shown, not assumed.

## 4. What is already known

From the round 3c data and the driver files:

| Item | Value | What it means for the network |
|---|---|---|
| Bass pair (two woofers in parallel) | 3.26 Ω minimum, at 147 Hz | Already below 3.5 Ω before the coax branch is added. The bass inductor's DCR helps a little. |
| Coax HF impedance | 24.7 Ω peak at 750 Hz (resonance), 8.7 Ω at 1.2 kHz | The high-pass sees a rising load near the resonance. |
| Sensitivity | coax LF 96 dB, bass pair about 90.5 dB, coax HF 112 dB | The coax LF needs about 5.5 dB of padding and the HF about 8–15 dB. |
| Active EQ the passive network can't copy | coax LF +10.6 dB at 1.83 kHz, bass delay 1.36 ms | The coax LF dip near 1.9 kHz stays partly unfilled. |
| SA550 drive | 25 V rms in the passive model (`[drive]` in the params) | Sets resistor power ratings and bass excursion. |

## 5. How it works

### 5.1 Simulation

- Each way's acoustic response is already known as pressure for 2.83 V at its terminals, in every
  direction (`speaker.system.Ways`). The passive network changes the voltage at those terminals. So
  each way's response is scaled by its network's voltage transfer, computed with the driver's complex
  impedance as the load (the `.zma` data from `speaker.vituixcad.impedances`).
- The networks are ladder circuits of R, L and C. Inductors carry their DCR and capacitors are ideal
  *(my call: capacitor ESR is too small to matter here)*.
- The whole speaker's input impedance is the bass branch and the coax branch in parallel, each with its
  network.
- The summed response then goes through the existing spinorama and `score_all`, unchanged.

### 5.2 Topology

The main spec's split stays: separate sub-networks joined at the terminals, with links that bypass the
parts the active modes don't use.

- **Bass low-pass:** series inductor and shunt capacitor (2nd order), with a 3rd-order option
  (an extra series inductor). An impedance-flattening RC (Zobel) across the woofers is optional.
- **Coax high-pass (bass-to-coax):** series capacitor and shunt inductor (2nd order). In semi-active
  mode a link shorts the capacitor and opens the shunt inductor, and the DSP does this crossover.
- **Coax LF low-pass:** 2nd or 3rd order, with an optional notch for the 2.7 kHz cone breakup.
- **Coax HF high-pass:** 3rd order (C, L, C), steep enough to keep drive out of the 750 Hz resonance,
  plus an L-pad for the 8–15 dB of attenuation and an optional RLC shunt to tame the 2–3.3 kHz hump.
- **Coax LF level:** a series resistor or L-pad for the ~5.5 dB match to the bass.
- The search tries each optional part on and off and keeps it only if it earns its place in the scores.

### 5.3 The fit

- **Target:** the same as the active fit, so the two are compared fairly: the summed response at 23°
  flat at the coax LF level, with the crossover points from the same grid (bass-to-coax 150–250 Hz;
  LF-to-HF 1200–2000 Hz, never below the driver's 1200 Hz).
- **Objective and selection** *(my call, as in the active fit)*: at each grid point pair the values
  are fitted for flatness at 23° alone. The other scores (window spread, presence lift, DI smoothness,
  whole-system rear rejection) only choose between grid point pairs, through `pick`.
- **Polarity** *(my call)*: each branch may be wired inverted. It's just how the wires go, so the
  search tries both for the bass and the coax HF.
- **Fit:** for each grid point pair, the component values are optimised continuously (least squares on
  the 23° response over the trusted band, with the minimum impedance as a hard constraint). Then each
  value snaps to a standard value and the rest are refitted around it: capacitors E12, resistors E24,
  inductors in 0.1 mH steps below 1 mH and 0.5 mH steps above *(my call)*.
- **Above the trusted limit,** only the HF's broad parts act (the L-pad and the 3rd-order slope). The
  optional RLC shunt is kept to a Q of 2 or less.
- **Picking the points:** the existing rule (`xo_fit.pick`): keep candidates within 1 dB of the best
  whole-system rear rejection, then take the first in directivity order.
- **Tolerance:** 200 random draws at ±5 % on every L and C. The report gives the spread of flat23,
  whole_rear and the minimum impedance. *(My call on pass and fail:)* the network fails if any draw
  takes the minimum impedance below 3.3 Ω, or if the 95th-percentile draw makes flat23 more than
  0.5 dB worse than the nominal values. A failing network goes back to the fit with the 3.5 Ω floor
  raised by 0.1 Ω at a time, up to 3.8 Ω.
- **The 3.5 Ω floor** applies to the modelled impedance at nominal values *(my call)*. The model is
  unmeasured, so the report also says how close the result comes to the floor, and the margin is
  checked again once the drivers are measured.

### 5.4 Semi-active mode

*(My calls.)*
- The coax boards are the same parts as in passive mode. In semi-active mode the links bypass the coax
  high-pass, and the DSP does the bass-to-coax crossover. The coax LF and HF networks and the coax pad
  are not refitted. One set of boards has to serve both modes.
- The DSP side is fitted with the existing active machinery, reduced to two channels:
  - an LR4 split on the same grid of bass-to-coax points;
  - bass gain, delay, polarity and up to 2 PEQs;
  - coax gain only, since its EQ is in the passive network.
  
  The DSP can delay the bass, so semi-active should recover most of the active rear rejection. The
  report shows how much.
- In semi-active mode the bass channel drives the bare woofer pair (3.26 Ω minimum) from its own amp.
  The report warns about this. It is not a pass or fail, because the amp isn't chosen yet.
- **The active comparison** is re-fitted in the same command with `xo_fit.fit_active`, from the same
  cached fields. It is deterministic and takes about 16 s, so the three modes are always compared on
  one set of data.

### 5.5 Parts

*(My calls, cost-aware; the owner can override.)*
- Air-core inductors in the coax network, where values are small and distortion matters more.
- Cored inductors (laminated or ferrite) for the bass low-pass, whose values (several mH) would be
  large, costly and lossy in air core. Their DCR is taken from typical catalogue values for the size.
- Polypropylene film capacitors in the coax network; bipolar electrolytics allowed above 30 µF, which
  in practice means the bass low-pass capacitor and any Zobel.
- Resistor power ratings for a sine at one third of the SA550's full voltage (about a tenth of full
  power, roughly what music averages), rounded up to the next standard rating. The full-power figure is
  listed beside it.

## 6. Outputs

In `out/<run>/passive/`, from a new command `speaker passive --params <file> --out <dir>`:

- `network.txt`: per mode, the netlist (part, value, type, tolerance, DCR, power rating), and which links
  to fit for each mode.
- `parts.csv`: the parts list for a pair.
- `semi_active.txt`: the DSP sheet for the two channels in semi-active mode, in the format of
  `crossover.txt`.
- `schematic_*.svg`: one schematic per sub-network, plus a link diagram for the three modes. It shows
  which parts each mode bypasses and the terminal pairs from main spec §7.2: one pair for passive, two
  for semi-active and three for fully active. Where the terminals go on the cabinet, and what kind
  they are, stays with the cabinet-construction work. *(My call: drawn with `schemdraw`, a new
  pure-Python dependency.)*
- `report.html`: scores for passive, semi-active and active side by side; 23° responses; polar maps;
  rear rejection; the impedance with the 3.5 Ω line; the tolerance spread; and where passive falls
  short of active, in plain sentences.
- `vituixcad_netlist.txt`: the network, laid out in the order you'd enter it in VituixCAD.

The owner's report page (the artifact) gets a "Passive" option next to "Adopted" and "Previous fit".
A round write-up goes in `docs/rounds/round-3d.md`.

## 7. What happens if the load check fails

If no network keeps the speaker at or above 3.5 Ω with the woofers in parallel, the report says so,
and gives the best network for wiring the woofers in series instead. In series, the pair is 12.8 Ω but
6 dB quieter (about 84.5 dB), so the coax must be padded about 6 dB more and the passive version plays
about 6 dB less loud. The owner then chooses. *(Nothing is decided here: it depends on the numbers.)*

## 8. Left to the builder

- The optimiser and its starting values.
- How the optional parts are searched (on/off combinations or a greedy add).
- Plot styling, within the existing report's look.
