# Round 3: coax offset against the mouth, top margin and bass heights

Date: 2026-10-28/29. The third optimisation round (spec §2), based on the
round 2 recommendation, n12-sq-m320 (`params/round-2.toml`). One shape
stage tried the coax offset against the flare mouth, a taller cabinet above
the coax, and higher bass drivers. The best two went to full-resolution
finals.

## What was run

- **Code commit:** `27b05eb` (branch `FyrBridd/round-3`).
- **Round file:** `rounds/round-3.toml`, sha256 `74fa4dd0…`.
- **BEM code version** (cache key): `d00d81a6…`, the same as rounds 1 and 2.
- **Command:** `uv run speaker round rounds/round-3.toml`, then the mesh
  check of `bass260-580`.
- **Outputs** (not committed): `out/rounds/3/report.html`, with
  `notes.md`, `results.json`, `finals.json`, each variant's files and the
  finals' own reports.

| Stage | Settings | Ranked by | Variants | Rejected | Errors | Time |
|---|---|---|---|---|---|---|
| shape | screen (5 per λ, to 2715 Hz) | directivity, behind the cardioid gate | 8 + stage base + round 2 | 0 | 0 | 3 h 23 min |
| finals | full run (6 per λ, with the puck) | – | 2 | – | 0 | 44 min |
| Total | | | | | | 4 h 07 min (22:22–02:29) |
| mesh check | bass260-580, 6 against 8 per λ | – | 1 | – | 0 | 1 h 17 min |

There was no vents stage. The rear rejection beat the 10 dB gate by about
7 dB, and moving the vents forward in round 2 cost 23° flatness. Every
variant was still tuned (felt resistance and chamber volume).

The stage base tuned to the round 2 design's own values (R 95.8 rayl,
V 16.0 L), so it is identical to the "round0" reference row, and its final
reused round 2's cached full run (11 s).

## The screen

Screen scores, trusted, in dB:

| Rank | Variant | Window spread | Presence lift | Sum | DI (bass + LF–HF) | flat23 | Hard | R / V |
|---|---|---|---|---|---|---|---|---|
| 1 | stage base (round 2 design) | 3.26 | 3.45 | 6.71 | 1.81 | 3.17 | 17.15 | 96 / 16.0 |
| 2 | bass260-580 | 3.27 | 3.45 | 6.71 | 1.82 | 3.01 | 17.15 | 96 / 16.0 |
| 3 | bass280-600 | 3.27 | 3.45 | 6.71 | 1.82 | 3.00 | 17.15 | 96 / 16.0 |
| 4 | in5-m310 | 3.19 | 3.44 | 6.63 | 1.85 | 3.15 | 17.15 | 96 / 16.0 |
| 5 | in10-m300 | 3.16 | 3.45 | 6.61 | 1.92 | 3.17 | 17.12 | 96 / 16.0 |
| 6 | in10-m290 | 3.18 | 3.46 | 6.64 | 1.98 | 3.20 | 17.12 | 96 / 16.0 |
| 7 | top60 | 3.47 | 3.50 | 6.97 | 1.40 | 3.33 | 17.05 | 115 / 13.5 |
| 8 | in15-m290 | 3.29 | 3.47 | 6.76 | 1.98 | 3.29 | 17.13 | 96 / 16.0 |
| 9 | top100 | 3.57 | 3.50 | 7.07 | 1.41 | 3.35 | 16.77 | 96 / 16.5 |

Ranking is window spread plus presence lift in 0.25 dB buckets, then DI
smoothness, then 23° flatness. Ranks 1–6 share a bucket (6.50–6.75 dB), so
the DI decides them.

What this shows:

- **Coax offset inwards is a tie.** 5–10 mm gains 0.07–0.10 dB of window
  spread, but loses 0.04–0.15 dB of LF–HF DI smoothness. With the mouth cut
  to fit (the flat front must hold the mouth + 2 × the offset + 2 × 5 mm),
  15 mm is worse on both. An offset would also make the pair handed.
- **A taller cabinet trades window for DI.** 60 or 100 mm above the coax
  (against 20 mm) costs 0.26–0.36 dB of window and sum, and gains 0.3 dB of
  LF–HF DI smoothness. The taller chamber needs more fillers: the tune
  picked 13.5 L at 60 mm (9.7 L of fillers) and 16.5 L at 100 mm (10.1 L).
- **Higher bass drivers help 23° flatness** by 0.16 dB on the screen and
  change nothing else by more than 0.01 dB. 280/600 mm is no better than
  260/580 mm.
- **The rejection does not move.** Every variant is within 0.4 dB of the
  base on hard rejection, and the felt-spread worst case stays at 13.5–13.7
  dB.

## Finals

In dB, trusted unless marked:

| Score | stage base (round 2 design) | bass260-580 |
|---|---|---|
| cardioid_hard (150–300) | 17.15 | 17.15 |
| cardioid_aim (300–500) | 12.35 | 12.36 |
| window_spread (150–2851) | 3.48 | 3.48 |
| presence_lift (2000–2851) | 4.10 | 4.10 |
| presence_lift (2851–5000, low) | 6.08 | 6.08 |
| flat23 (150–2851) | 4.21 | **4.01** |
| flat23 (above, low) | 10.45 | 10.45 |
| di_smooth_bass | 0.57 | **0.55** |
| di_smooth_lfhf | 1.26 | 1.26 |
| whole_rear (higher is better) | 3.45 | **3.70** |

Crossover gains: bass +11.5 dB (round 2: +11.7), coax LF 0, coax HF −4.8 dB.

bass260-580 is equal or better on every row. The low-confidence rows above
the trusted limit are identical, so it has no above-limit penalty like round
2's cone.

## Recommendation

**Owner's decision (2026-10-29): bass260-580** (`params/round-3.toml`).
The first record of this decision, dcd043a, said the owner kept
n12-sq-m320; that was a mistyped design name.

**bass260-580**: the round 2 recommendation with the bass driver centres at
260 mm and 580 mm (were 200 mm and 520 mm). Everything else is unchanged.

This is not the round report's recommendation. The report takes the first
final, the stage base, which won the screen's DI tie-break by 0.01 dB. On
the finals bass260-580 wins that tie-break by 0.02 dB and is ahead on 23°
flatness. The finals are not re-ranked (M4); this is the second round
running where that picked the wrong final.

The gain is small but free: 0.2 dB of 23° flatness, 0.25 dB more whole-
speaker rear rejection at 150–300 Hz, the same bass volume, and 0.17 dB
less bass gain.

## Trusted limit of the recommendation

**2851 Hz**, the same as round 2. The mesh check of bass260-580 (6 against
8 elements per wavelength, rejection capped at 30 dB) agrees within 0.06 dB
in every band up to 2851 Hz (worst in the 1 kHz band). It stopped at the
4 kHz band on memory (about 25.2 GB against 22 GB).

## Still open

1. **Shape search has converged.** Three rounds of front and cabinet
   variants now move the window by about 0.1 dB at a time. The 2–2.85 kHz
   window rise (presence lift 4.1 dB, warn) has not responded to any of
   them. It needs the crossover and EQ (plan 3).
2. **Re-rank the finals on their own scores (M4).** It picked the wrong
   final in rounds 2 and 3.
3. **The drivers are not bought.** The provisional assumptions of round 0
   (cone depth 65 mm, HF exit 25.4 mm, cutout 230 mm, fixed crossover,
   basket-cylinder coax) still stand, and the fillers will need re-sizing
   once the coax is measured.
4. **Round 0's open questions 3–5** (the rear null at ±120°, the low bass
   below the axis, the puck detail above the trusted limit) are unchanged.
