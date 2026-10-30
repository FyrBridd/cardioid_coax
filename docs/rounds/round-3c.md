# Round 3c: the crossover search re-ranked

Date: 2026-10-29. The owner asked why the fit chose a 1.0 kHz LF-to-HF point when BMS recommends
1.2 kHz, and approved four changes to the ranking. This round reruns the design and round 3 under
them. Everything here is simulated; the drivers are not bought.

## What changed (scoring version 3)

1. **DI smoothness over fixed bands.** di_smooth_bass is scored over 120–500 Hz and di_smooth_lfhf over
   500 Hz–4 kHz (trusted to 2851 Hz) for every crossover. Before, each candidate was scored over fc/2 to 2fc of its
   own point. So 1.0 kHz was judged on 500–2000 Hz, and 1.2 kHz on 600–2400 Hz, which takes in the DI step above
   2 kHz. The bass band starts at 120 Hz, where the screen grid starts, so screens and finals score it
   alike.
2. **Whole-system rear rejection in the crossover search.** Candidates more than 1 dB below the best
   trusted whole_rear (150–300 Hz) are dropped; the rest are ranked as before. The search had ignored
   it, though it is the point of the cardioid.
3. **The driver's floor.** The LF-HF grid stops at the coax datasheet's recommended HF crossover,
   `hf_xo_min_hz = 1200` in `data/drivers/bms_10c262.toml` (BMS 10C262 datasheet, HF section).
4. **Tiers, not fixed buckets.** A tie on window_spread + presence_lift is any value within 0.25 dB of
   the best one left. The old buckets were fixed multiples of 0.25 dB and split values 0.02 dB apart.
   This applies to the round ranking and the crossover search.
5. **Finals on their own scores (M4).** The round report ranks the finals on their own finer-mesh scores and
   recommends the best one, not the first in screen order. A final without its own scores, or one that fails the
   gate, is tagged out.

The spec is amended to match (`docs/superpowers/specs/2026-10-26-directivity-optimisation-design.md`,
§4 and §5–6). Stored round results are re-scored because the scoring version changed.

## What was run

- **Code:** branch `FyrBridd/xo-ranking`, commit `f70d5e6`, scoring version 3. The BEM code and cache
  are unchanged, so there were no new solves. Outputs are not committed.
- **Full run:** `uv run speaker directivity --params params/round-3.toml --out out/xo-ranking/round-3-design --trusted-limit 2851`,
  16 s. It wrote `crossover.txt`, `crossover.json` and the VituixCAD export (`vituixcad/`).
- **Round:** `uv run speaker round rounds/round-3.toml --out out/rounds/3c`, 3 min 35 s. All 10 shape variants were
  re-scored, and both finals took 14 s each from the cache.
- **Tests:** `uv run pytest -q` gives 372 passed, 13 deselected.

## The crossover of bass260-580 (`params/round-3.toml`)

The search now picks **bass to coax 150 Hz, coax LF to HF 1200 Hz**. The filters are identical to the
fit with the grid floored at 1200 Hz that the owner was shown before approving.

| Way | Gain | Delay | Polarity | PEQ (f0, gain, Q) |
|---|---|---|---|---|
| bass | +11.5 dB | 1.36 ms | inverted | 80 Hz +1.2 dB Q1.6; 214 Hz −1.8 dB Q1.2 |
| coax LF | 0 dB | 0.01 ms | normal | 153 Hz +8.9 dB Q1.8; 544 Hz −4.1 dB Q0.9; 1425 Hz −3.4 dB Q3.8; 1700 Hz +12.0 dB Q0.8 (limit) |
| coax HF | −8.0 dB | 0 | normal | 1635 Hz +3.0 dB Q1.0; 1749 Hz −4.2 dB Q1.8; *low:* 7505 Hz +9.0 dB Q3.9; 10842 Hz +12.0 dB Q3.2 (limit); 14096 Hz +12.0 dB Q3.0 (limit) |

Against the round 3b fit (230 / 1000 Hz), both scored here with version 3. The 3b fit is rebuilt with the same
filter fit, and it reproduces the 3b scores.

| Trusted score (dB) | 150 / 1200 Hz | 230 / 1000 Hz |
|---|---|---|
| Whole-system rear rejection, 150–300 Hz (higher is better) | **9.88** | 4.16 |
| 23° flatness | 1.16 | **0.92** |
| Window spread | **3.40** | 3.45 |
| Presence lift, 2–2.85 kHz | 3.90 | **3.64** |
| DI smoothness, bass, 120–500 Hz | 1.66 | **1.32** |
| DI smoothness, LF–HF, 500–2851 Hz | **1.08** | 1.19 |
| Cardioid hard / aim (coax section) | 17.15 / 12.36 | 17.15 / 12.36 |

The 150 Hz point buys 5.7 dB of whole-system rear rejection. It costs 0.24 dB of 23° flatness, 0.26 dB of
presence lift and 0.34 dB of bass DI smoothness. Window plus presence is 7.30 against 7.09 dB, within one
0.25 dB tier. Part of the presence gap is the low-confidence HF filters: on the trusted-band data the
search ranks on, it is 7.13 against 7.16 dB.

The raw coax LF and HF sound-power DI agree within 0.4 dB from 700 Hz to 1.3 kHz, because the horn fires
through the cone. So the LF-HF point itself barely changes the directivity. The 1.2 kHz point gives the
compression driver 6 dB less signal at its 750 Hz resonance (−25 against −19 dB) and 4.5 dB less at 1 kHz.

**Excursion** (lumped, sine, passband level at 2 m): the woofers still set the limit, at 103.0 dB at 20 Hz,
105.5 dB at 40 Hz, and 107.5 dB at 40 Hz with a 35 Hz 12 dB/oct subsonic high-pass. The coax LF's +8.9 dB at
153 Hz brings its xmax point down to 117.7 dB at 154 Hz (from 128.7 dB), still about 15 dB clear of the
woofers.

**DSP sheet:** use 96 kHz. At 48 kHz the biquads miss the coax HF by 0.50 dB at 2–5 kHz and the coax LF by
0.25 dB. The largest gain plus PEQ boost is +12.5 dB, the bass at 79 Hz.

**Passive-network warnings:**
- **EQ boost:** the coax LF needs 10.6 dB at 1.83 kHz and the coax HF 16.8 dB at 11 kHz (passive: 6 dB).
- **HF resonance:** 1200 Hz is within an octave of the 750 Hz HF resonance. The warning's limit is 1500 Hz;
  BMS recommends 1200 Hz.
- **Low-confidence filters:** the three coax HF filters above 2851 Hz rest on low-confidence data. As in
  round 3b, leave them out until the coax is measured.

## Round 3c

The shape stage re-ranked under version 3: bass280-600, bass260-580, stage-base, round0, in5-m310,
in10-m290, top60, in10-m300, top100, in15-m290. None failed the gate. Screen results do not record their
crossover points; the top four share their window, presence and flatness scores. The bass heights move only
the bass DI and whole_rear.

The finals, ranked on their own scores (both at 150 / 1200 Hz):

| Final | Window + presence | Bass DI | LF–HF DI | flat23 | whole_rear |
|---|---|---|---|---|---|
| bass280-600 | 7.31 | 1.60 | 1.08 | 1.17 | 9.92 |
| bass260-580 | 7.30 | 1.66 | 1.08 | 1.16 | 9.88 |

The report picks bass280-600 (bass centres 280 / 600 mm, 20 mm higher) on 0.06 dB of bass DI. That is well
inside what the model can resolve. **I recommend keeping bass260-580**, the owner's choice, so the design file,
the CAD and the exports stay as they are. The round's mechanical pick is noted, not adopted.

## Open

- The owner hasn't checked the VituixCAD import yet. The new export is in
  `out/xo-ranking/round-3-design/vituixcad/`.
- The fully active setup needs six amplifier channels per pair: bass, coax LF and coax HF for each speaker.
- The coax HF filters above the trusted limit wait for a measurement.
