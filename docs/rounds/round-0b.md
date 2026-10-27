# Round 0b: blend moved above the trusted limit

Date: 2026-10-27. Same design as round 0. Only the way the cabinet and puck
models are joined has changed. Nothing was optimised.

## What changed

Your choice on round 0's observation 3 (2026-10-27): the cabinet-to-puck
crossfade now starts at the trusted limit instead of ending there
(commits `dd5300f`, `1e5e464`, branch `FyrBridd/blend-above-limit`).

- **Before:** the fade ran 1425–2851 Hz, inside the trusted band.
- **Now:** it runs **2851–3805 Hz** (0.42 octave), above the trusted limit,
  up to the top of the cabinet's 4 kHz band. The vents fade out across it.
- So every trusted score is now the checked cabinet BEM alone. The cabinet
  fields inside the blend (above 2851 Hz) are not mesh-checked, and the
  report says so.
- The 4 kHz cabinet band, unused in round 0, is now used by the blend.

## What was run

- **Code:** `1e5e464`. Params file unchanged (sha256 `f877a7c7…`).
- **BEM code version** (cache key): `daf4c4ca…`. It differs from round 0's
  `f6c88fea…` because gmsh is now 4.15.2 from PyPI, so every band was
  re-solved.
- **Commands:** as round 0 (mesh check, then the full run).

## Trusted 3D limit: still 2.85 kHz

The mesh check on the new gmsh gives the same result: every checked band
agrees within 0.08 dB (worst: the 160 Hz band), and the 4 kHz check mesh
still needs about 24 GB against the 22 GB budget.

## Scores that changed

Only scores whose band touched the old fade moved. The cardioid, bass and
rear rows are identical to round 0.

| Key | Band (Hz) | Confidence | Round 0 | Round 0b | Puck share then → now |
|---|---|---|---|---|---|
| presence_lift | 2000–2851 | trusted | 5.67 dB warn | **5.26 dB warn** | 100% / 0.72 → 0% |
| window_spread | 150–2851 | trusted | 4.70 dB warn | **4.05 dB warn** | 24% / 0.13 → 0% |
| flat23 | 150–2851 | trusted | 4.10 dB | 4.14 dB | 24% / 0.13 → 0% |
| di_smooth_lfhf | 700–2800 | trusted | 1.42 dB | 1.60 dB (SPDI 1.60, ERDI 0.69) | 48% / 0.24 → 0% |
| rolloff_monotonic | 150–2851 | trusted | 0 pass | 0 pass | 24% / 0.13 → 0% |
| flat23 | 2851–16000 | low | 9.57 dB | 9.55 dB | 100% → 97% / 0.90 |
| presence_lift | 2851–5000 | low | 7.15 dB | 7.06 dB | 100% → 91% |
| window_spread | 2851–16000 | low | 15.37 dB | 15.37 dB | 100% → 97% / 0.90 |

What this means:

- **The puck was making the presence region look worse.** With the checked
  cabinet model, the trusted presence lift drops by 0.4 dB and the window
  spread by 0.65 dB. Both still warn. The window stays within ±1 dB only up
  to 898 Hz, as before; above that the 0–10° directions rise to about +3.5–4 dB
  over 23° by 2–2.8 kHz.
- **The DI at the LF–HF crossover is a little less smooth** (1.42 → 1.60 dB).
  The puck's round edges had smoothed it slightly.
- **The low-confidence scores barely move.** Above 3.8 kHz they are puck
  only, as before.

## Plots

- The hatched blend now sits at 2851–3805 Hz on the spinorama, polar maps
  and window family.
- In the blend the on-axis and 0–10° curves still show a peak near 3.3 kHz
  and a dip of about 5 dB near 3.6 kHz. Round 0 showed the same shape from
  the puck alone (its observation 5), so it is not a crossfade artefact.
  It is still low confidence.
- 1.4–2.85 kHz is now cabinet only. The shape is close to round 0's: on-axis
  peaks near 2 kHz at about 101 dB at 2 m.

## Still open

Round 0's open questions 1–5 stand: the cardioid fails its hard gate
(5.3 dB), the vents add to the rear sound above about 180 Hz, the coax LF
peak at 550–650 Hz, the 1.8–1.9 kHz notch, and the puck's exaggerated
detail above the trusted limit. Question 6 (the 2 m distance) is settled:
you confirmed 2 m is fine.

## Run times

| Step | Start | End | Time |
|---|---|---|---|
| Mesh check | 06:11 | 07:49 | 1 h 38 min |
| Full run | 07:49 | 09:20 | 1 h 30 min |
| Total | | | 3 h 08 min |

Element counts are as in round 0. The 12 kHz puck band took 3242 s.
