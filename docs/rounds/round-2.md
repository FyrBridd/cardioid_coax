# Round 2: vents forward, then the flare and the front

Date: 2026-10-28. The second optimisation round (spec §2), based on the round
1 recommendation (`params/round-1.toml`). Stage 1 moved the vents further
forward on shorter slots, tuning the felt and the coax chamber volume for
each. Stage 2 tried the flare mouth, small edges with the full mouth, a
deeper flare, and the six rows round 1 dropped, on the stage 1 winner. The
best two went to full-resolution finals. One more combination was run
separately afterwards (below).

## What was run

- **Code commit:** `74da03c` (branch `FyrBridd/round-2`).
- **Round file:** `rounds/round-2.toml`, sha256 `8cac36bc…`.
- **BEM code version** (cache key): `d00d81a6…`, the same as round 1.
- **Command:** `uv run speaker round rounds/round-2.toml`, then a full run
  and mesh check of `n12-sq-m320` (below).
- **Outputs** (not committed): `out/rounds/2/report.html`, with
  `notes.md`, `results.json`, `finals.json`, each variant's files, the
  finals' own reports, and `extra/n12-sq-m320/`.

| Stage | Settings | Ranked by | Variants | Rejected | Errors | Time |
|---|---|---|---|---|---|---|
| vents | vents (to 600 Hz) | tuned cardioid margin | 7 + stage base + round 1 | 0 | 0 | 51 min |
| shape | screen (5 per λ, to 2715 Hz) | directivity, behind the cardioid gate | 13 + stage base + round 1 | 1 | 0 | 4 h 08 min |
| finals | full run (6 per λ, with the puck) | – | 2 | – | 0 | 4 h 06 min |
| Total | | | | | | 9 h 05 min (10:32–19:37) |
| extra | full run, then mesh check | – | 1 | – | 0 | 45 min + 1 h 35 min |

The report calls the base design "the round 1 recommendation" (the new
`base_label` round-file key, da30def). The "round0" rows in each stage are
that design, untuned.

Two changes from the proposal: the short slots sit at 88 mm, not 85 mm,
because the fillers need the slots to start at least 51 mm behind the front;
and the 54 mm outer baffle was replaced by 30 mm, because 54 mm leaves the
bass section at 48.2 L, under its 50 L minimum.

## Vents further forward

The 100 mm slots cannot start closer than 51 mm behind the front, so they
stop at 105 mm. The 70 mm slots go to 88 mm.

| Variant | Slots (from front, width, per side) | R (rayl) | V (L) | Hard | Aim | Worst-case margin | flat23 (screen) |
|---|---|---|---|---|---|---|---|
| s70-cff=88-n=9 (winner) | 88 mm, 70 mm, 9 | 115 | 10.0 | 17.66 | 12.29 | 3.70 | 3.80 |
| s70-cff=88-n=12 | 88 mm, 70 mm, 12 | 96 | 16.0 | 17.09 | 12.28 | 3.66 | 3.48 |
| s70-cff=95-n=12 | 95 mm, 70 mm, 12 | 96 | 16.5 | 17.03 | 12.13 | 3.56 | 3.34 |
| s70-cff=95-n=9 | 95 mm, 70 mm, 9 | 96 | 12.5 | 16.76 | 12.10 | 3.44 | 3.50 |
| 105 mm, 9 slots | 105 mm, 100 mm, 9 | 115 | 15.5 | 16.80 | 11.87 | 3.28 | 3.25 |
| 105 mm, 12 slots | 105 mm, 100 mm, 12 | 137 | 17.0 | 16.71 | 11.82 | 3.23 | 3.28 |
| 105 mm, 6 slots | 105 mm, 100 mm, 6 | 115 | 10.5 | 16.46 | 11.82 | 3.09 | 3.47 |
| stage base (round 1) | 120 mm, 100 mm, 9 | 115 | 16.5 | 16.27 | 11.56 | 2.95 | 3.00 |

Hard is 150–300 Hz, aim 300–500 Hz, both at the tuned R and V. The margin
is the worst over the ±30 % felt spread, above the 10 / 6 dB gates.

- **Forward still helps, but less.** From 120 mm to 88 mm the best margin
  rises 0.75 dB, against about 1.5 dB per 40 mm in round 1.
- **The two 88 mm layouts tie.** 9 slots won by 0.04 dB of margin, but its
  tune cut the coax chamber to 10 L. 12 slots keeps 16 L.
- **The forward vents cost 23° flatness.** Every forward layout is 0.25–0.8
  dB less flat at 23° than the stage base on the vents screen; the 10 L
  chambers are the worst.
- **Felt spread** (12-slot, 88 mm): hard / aim = 13.8 / 13.0 dB at R × 0.7,
  13.7 / 11.7 dB at R × 1.3 (round 1: 13.5 / 12.2 and 12.9 / 11.0).

## The flare and the front

All shape variants carry the vents winner (9 slots, 10 L), except the two
deeper flares, which need the 100 mm slots at 120 mm. Screen scores,
trusted, in dB:

| Rank | Variant | Window spread | Presence lift | Sum | DI (bass + LF–HF) | flat23 | Hard |
|---|---|---|---|---|---|---|---|
| 1 | cone (cone flare) | 3.24 | 3.38 | 6.61 | 1.31 | 3.60 | 17.81 |
| 2 | sq-m320 (320 mm mouth) | 3.27 | 3.45 | 6.72 | 1.67 | 3.62 | 17.75 |
| 3 | in10 (coax 10 mm inwards) | 3.19 | 3.47 | 6.66 | 1.83 | 3.60 | 17.70 |
| 4 | sq-m310 | 3.30 | 3.49 | 6.78 | 1.71 | 3.64 | 17.72 |
| 5 | stage base (290 mm mouth) | 3.32 | 3.53 | 6.86 | 1.83 | 3.66 | 17.66 |
| 6 | bass260-580 (lower bass) | 3.32 | 3.53 | 6.86 | 1.85 | 3.46 | 17.66 |
| 7 | sq-m270 | 3.32 | 3.53 | 6.85 | 2.01 | 3.68 | 17.59 |
| 8 | round 1 (reference) | 3.29 | 3.50 | 6.79 | 2.07 | 3.05 | 16.27 |
| 9 | ro15-m290 (15 mm roundover) | 3.55 | 3.66 | 7.22 | 1.86 | 3.70 | 17.59 |
| 10 | out10 (coax 10 mm outwards) | 3.61 | 3.66 | 7.27 | 1.87 | 3.73 | 17.70 |
| 11 | ch15-m290 (15 mm chamfer) | 3.73 | 3.78 | 7.51 | 1.88 | 3.73 | 17.50 |
| 12 | ob30 (30 mm outer baffle) | 4.30 | 4.27 | 8.57 | 2.48 | 3.62 | 15.81 |
| 13 | ob36 (36 mm outer baffle) | 4.80 | 4.70 | 9.50 | 2.64 | 3.96 | 15.76 |
| – | top60 (60 mm top margin) | rejected | | | | | |

Ranking is window spread plus presence lift in 0.25 dB buckets, then DI
smoothness, then 23° flatness; that is why in10 ranks below sq-m320.

What this shows:

- **The edge itself hurts.** A 15 mm roundover or chamfer with the full
  290 mm mouth is 0.35–0.65 dB worse than the square edge. So round 1's edge
  result was the edge, not the smaller mouths it was paired with.
- **A wider mouth helps a little.** 270 / 290 / 310 / 320 mm give 6.85 /
  6.86 / 6.78 / 6.72 dB.
- **A deeper flare hurts a lot.** The 30 and 36 mm outer baffles are 1.7–2.6
  dB worse on the window. Their hard rejection (15.8 dB) is also 0.5 dB
  below round 1's 16.3 dB on the same 120 mm vents.
- **The coax offset matters.** 10 mm inwards gains 0.2 dB; 10 mm outwards
  loses 0.4 dB.
- **The lower bass heights** change only the 23° flatness (0.2 dB better)
  and the DI by 0.02 dB.
- **top60 was not tested.** It was rejected before tuning: "rear block would
  need to be 357.9 mm tall to fill 13.10 L, but the coax chamber is 330.0 mm
  tall; raise the target net volume". The carried 10 L does not fit the
  shorter chamber, and the layout check runs before the tune could pick a
  volume that does (deferred minor M2).
- **The new vents did not change directivity.** Round 1 untuned and the
  stage base differ by at most 0.07 dB on window spread and presence lift,
  but round 1 is 0.6 dB flatter at 23°.

## Finals

The finals rerun the top two at 6 elements per wavelength, with the puck and
the full 16 kHz grid. In dB, trusted unless marked:

| Score | cone | sq-m320 | Round 1 final |
|---|---|---|---|
| cardioid_hard (150–300) | 17.81 | 17.75 | 16.27 |
| cardioid_aim (300–500) | 12.47 | 12.37 | 11.56 |
| window_spread (150–2851) | 3.36 | 3.49 | 3.55 |
| window_spread (above, low) | 16.65 | 14.79 | 14.76 |
| presence_lift (2000–2851) | 3.96 | 4.11 | 4.20 |
| presence_lift (2851–5000, low) | **8.83** | 6.10 | 5.78 |
| flat23 (150–2851) | 4.52 | 4.62 | 3.60 |
| flat23 (above, low) | 13.16 | 10.45 | 9.79 |
| di_smooth_bass | 0.51 | 0.51 | 0.62 |
| di_smooth_lfhf | 0.87 | 1.26 | 1.44 |
| whole_rear | 3.31 | 3.30 | 3.58 |

The screen order held: the cone wins every trusted directivity score. But
it is clearly worse just above the trusted limit: 2.7 dB more presence lift
at 2.85–5 kHz, and 2.7 dB less flat at 23° there. Those rows are low
confidence and do not count in the ranking, but the difference is large and
plausible. The cone meets the baffle at a sharp corner, the ellipse does
not, and the puck model above the limit includes the flare profile
(`coax_geometry.py`). 2–5 kHz is also the band the owner is most sensitive
to.

Both finals carry the 10 L chamber and are 0.9–1.0 dB less flat at 23° than
round 1.

## The extra run: n12-sq-m320

Since the cone is ruled out on its 2.85–5 kHz rows, and the 10 L chamber
costs flatness and filler print, I ran the 320 mm mouth on the tied 12-slot
vents: `s70-cff=88-n=12`'s tuned params (R 95.8 rayl, V 16.0 L) with
`flare_mouth_d_mm = 320`. Full run (6 per λ, with the puck), then the mesh
check.

| Score | n12-sq-m320 | sq-m320 (9 slots) | Round 1 final |
|---|---|---|---|
| cardioid_hard | 17.15 | 17.75 | 16.27 |
| cardioid_aim | 12.35 | 12.37 | 11.56 |
| window_spread | 3.48 | 3.49 | 3.55 |
| presence_lift | 4.10 | 4.11 | 4.20 |
| presence_lift (2851–5000, low) | 6.08 | 6.10 | 5.78 |
| flat23 | 4.21 | 4.62 | 3.60 |
| di_smooth_bass | 0.57 | 0.51 | 0.62 |
| di_smooth_lfhf | 1.26 | 1.26 | 1.44 |
| whole_rear | 3.45 | 3.30 | 3.58 |

- 12 slots win back 0.4 dB of the 23° flatness for 0.6 dB of hard
  rejection, still 0.9 dB above round 1. Everything else is within 0.06 dB
  of the 9-slot version.
- **Fillers:** 3.7 L (16.0 L net from a 19.7 L empty chamber), against
  9.7 L for the 10 L designs and 3.2 L in round 1.
- **Crossover gains:** bass +11.7 dB, coax LF 0, coax HF −4.8 dB (round 1:
  +11.6, 0, −4.9).
- **Trusted limit: 2851 Hz.** The mesh check (6 against 8 per λ, rejection
  capped at 30 dB) agrees within 0.06 dB in every band up to 2851 Hz (worst
  at 2520 Hz). It stopped at the 4 kHz band on memory (about 25.2 GB against
  22 GB).

## Recommendation

**n12-sq-m320.** It is the round 1 recommendation with 70 mm slots (was
100 mm) centred 88 mm from the front (was 120 mm), 12 slots per side (was
9), R ≈ 96 rayl (was 115), V = 16.0 L (was 16.5 L), and a 320 mm elliptical
flare mouth (was 290 mm). Square edges, 330 mm wide, everything else
unchanged.

This is not the round report's recommendation. The report recommends the
first final, the cone, and shows its dimensions and curves. The finals are
not re-ranked, and nothing in the ranking looks above the trusted limit
(M4). The extra run's own report is `out/rounds/2/extra/n12-sq-m320/report.html`.

Against round 1 it is a modest step: 0.9 dB more hard rejection, 0.8 dB more
aim, a 13.7 dB felt-spread worst case against 12.9 dB, and smoother DI. The
window and presence lift are 0.1 dB better, which is not audible. It costs
0.6 dB of 23° flatness, from the forward vents. Round 1 already passed the
gate by 6.3 dB, so this buys margin rather than fixing a failure. Keeping
round 1 is a reasonable choice.

## Still open

1. **top60** (60 mm top margin) is untested after two rounds. Fix M2: let
   the tune choose a volume before the layout check rejects the carried one.
2. **Re-rank the finals** on their own scores and flag a reversal (M4), and
   give the 2.85–5 kHz rows some say in the ranking, even at low weight.
3. **Mouth beyond 320 mm, or a smoother flare profile.** The mouth gain is
   steady but small. Check how large a mouth the 330 mm front allows.
4. **in10 with the 320 mm mouth.** The inward offset gained about as much as
   the mouth; the two were not combined.
5. **The 2–2.85 kHz window rise** is still the main directivity problem. No
   front shape in two rounds removed it. It needs the crossover (plan 3).
6. **Round 0's open questions 3–5** (the rear null at ±120°, the low bass
   below the axis, the puck detail above the trusted limit) are unchanged.
   The provisional assumptions of round 0 (cone depth 65 mm, HF exit
   25.4 mm, cutout 230 mm, fixed crossover, basket-cylinder coax) still
   stand, and the drivers are not bought yet.
