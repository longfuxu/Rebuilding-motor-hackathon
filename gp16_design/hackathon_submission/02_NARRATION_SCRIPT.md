# Demo Video — Narration Script (Eleven Labs → screen-recorded deck)

**Runtime:** full version ≈ 3:00 (14 slides) · **2-minute cut** = drop the three `[CUT-FIRST]` blocks (slides 8, 11, 12).
**Pairs with:** `deck.html` (present full-screen; advance with → / Space).

> **✅ The audio is already generated** (high-quality neural voice — Microsoft *Andrew*, via `edge-tts`, free, no key needed). Files are in **`audio/`**: per-slide `v01.mp3 … v14.mp3`, plus three ready-to-drop combined tracks: `_full_all14.mp3` (4:56), `_cut_11slides.mp3` (4:01), `_cut_spine8.mp3` (**3:00 — recommended demo length**).
> - Regenerate / change voice: `python generate_voiceover_edge.py` (try `VOICE="en-US-BrianMultilingualNeural"` for a deeper voice, or `RATE="+8%"` to speed up).
> - **Eleven Labs** (`generate_voiceover.py`, best-quality config) is optional and currently **out of credits** on your key (5 left; ~55–102 needed per line). Top up or swap the key, then `python generate_voiceover.py` for Eleven Labs quality instead.

### How to use this file
1. In **Eleven Labs**, pick a calm documentary/science voice (*Adam*, *Daniel*, or *Charlie*). Settings: **Stability ~50, Similarity ~75, Style ~5–15, Speed 0.95–1.0.**
2. Paste **each numbered block's spoken text** and export one MP3 per slide — filenames below (`v01.mp3` … `v14.mp3`). Doing it per-slide makes syncing to slides trivial.
3. The spoken text is already written **phonetically** for TTS (e.g. "phi-29", "E-one-nineteen-Q"). Don't paste the on-screen cues or word counts. Pronunciation cheat-sheet at the bottom if a word comes out wrong.
4. Record: open `deck.html` full-screen, screen-record (QuickTime / OBS / Loom), advance each slide as its MP3 plays — or drop slide screenshots + MP3s on an iMovie/CapCut/Descript timeline, one per beat.
5. **Delivery:** unhurried, confident, land the *last sentence* of each block. Don't rush the numbers. Let the honesty beats (8 and 14) breathe.

---

## v01 — TITLE  ·  ~0:14
> Nature builds its strongest molecular motors as rings of identical parts. We rebuilt one as a single, addressable chain — to ask a question the ring itself hides. I'm Longfu Xu, from the Bustamante Lab at Berkeley.

## v02 — THE PROBLEM  ·  ~0:22
> The phi-29 gp16 motor packs DNA at nearly sixty piconewtons of force — as five identical subunits. But its real questions are per-subunit: which seat is the special, regulatory one? Is ATP fired one seat at a time, or all at once? In a ring of identical parts you can't tell. Put in a defect, and it lands at a random, unknown seat.

## v03 — THE MOVE  ·  ~0:18
> So we make the seats addressable. Fuse the five subunits into one gene, and a mutation in copy k lands at ring seat k — every time. It's the trick that made single-chain ClpX interpretable, twenty years ago. No one had done it for gp16. Rebuild it, to understand it.

## v04 — THE LEAD IS REAL  ·  ~0:19
> Here's the lead: cp233 — one seventeen-hundred-and-fifty-residue chain. Three independent structure predictors — Boltz-2, OpenFold3, and AlphaFold3 — all agree: five interfaces out of five. And it superimposes on the real motor — under two angstroms per subunit, all five copies landing on the native ring. The only extra density is the linker.

## v05 — PHYSICS  ·  ~0:17
> That's predictors. Now physics. In molecular dynamics, the apo ring stays shut, the ATP-bound native pops two seams open — and our design stays closed, all five catalytic interfaces engaged, with an interface energy even tighter than native. A physics engine that never saw our metrics, agreeing anyway.

## v06 — FUNCTION FILTER  ·  ~0:19
> But closing the ring is not enough. Drop each design into the real motor, with its DNA — and cp233's channel threads the double helix, just like native. Two rival cuts also close the ring, yet their channels are too narrow to thread DNA at all. Closing the ring is not the same as working as a motor.

## v07 — THE PAYOFF  ·  ~0:21
> And here's the payoff. We install a catalytically dead subunit — E-one-nineteen-Q — at one chosen seat, or at all five. The ring stays fully closed every time. A position control seals it: remove a single arginine finger, and it disables exactly the one interface we edited — and no other. The defect goes where we chose. That is genetic addressability.

## v08 — AN HONEST NULL  ·  ~0:16  `[CUT-FIRST]`
> One honest caveat. On gp16, the single chain is genuinely more coordinated than native — nearly one-point-eight fold. But test the same idea on ClpX, and that advantage vanishes. So this is not a universal law of covalency — it's a property of gp16's topology. We report the null as loudly as the win.

## v09 — IT GENERALIZES  ·  ~0:17
> Zoom out. Across thirteen ring motors, the whole build problem collapses to one rule: if a terminus fouls the substrate channel, use circular permutation; otherwise, fuse directly. Three for three on the systems where we already know the answer. From structure alone, we can predict how to build a motor into a single chain.

## v10 — VALIDATED ON A KNOWN ANSWER  ·  ~0:13
> And we tested the framework where the answer is already known. Run blind on ClpX, it passes the good motor six out of six — and catches the two dead-coupler mutants at zero out of six. It separates functional from broken.

## v11 — METHOD  ·  ~0:14  `[CUT-FIRST]`
> How? A five-level metric hierarchy that grades geometry and function — never global confidence. Plus an MSA-tiling trick that roughly doubled prediction confidence, and made a seventeen-hundred-residue ring foldable at all.

## v12 — WHAT IT UNLOCKS  ·  ~0:14  `[CUT-FIRST]`
> What it unlocks: set the ATP pattern at chosen seats, read out the force in optical tweezers, and finally decide whether firing is sequential or concerted. That's the bridge from design to measurement.

## v13 — HOW WE USED CLAUDE  ·  ~0:20
> All of this ran on Claude. Claude Code orchestrated every GPU job — and a workflow of nineteen sub-agents that scored, then adversarially re-checked, every fold. That verification changed the science: it caught false negatives, and it's why the honest nulls are in this talk at all. Claude Science ran the folding campaigns. Total compute: under ten dollars.

## v14 — CLOSE  ·  ~0:16
> Predict the build. Rank the coupling. Design something stronger — then measure it. It's cross-checked, not yet validated in the lab, and every prediction here is a falsifiable hypothesis. A molecular motor we can both predict and measure. Thank you.

---

## Timing summary (measured, neural voice)
| Version | Slides | Track file | Runtime |
|---|---|---|---|
| **Spine (recommended)** | v01,02,03,04,07,09,13,14 (8) | `audio/_cut_spine8.mp3` | **3:00** |
| 11-slide cut | drop v08, v11, v12 | `audio/_cut_11slides.mp3` | 4:01 |
| Full director's cut | all 14 | `audio/_full_all14.mp3` | 4:56 |

Per-slide durations (s): v01 16.6 · v02 26.6 · v03 19.9 · v04 27.0 · v05 22.8 · v06 23.0 · v07 26.1 · v08 24.0 · v09 22.3 · v10 14.7 · v11 17.1 · v12 13.6 · v13 24.5 · v14 17.2.

**Pick by your event's cap:** most Cerebral Valley demos are 2–3 min → use the **3:00 spine** (present slides 1,2,3,4,7,9,13,14 only, or hide the others). If you have up to 5 min, use the full track and all 14 slides. To fit all 14 into ~4 min, regenerate faster: `RATE="+15%" python generate_voiceover_edge.py`.

## Pronunciation cheat-sheet (if Eleven Labs mangles a term)
| Written | Say / type into Eleven Labs |
|---|---|
| φ29 / phi-29 | "fye twenty-nine" |
| gp16 | "gee-pee sixteen" |
| cp233 | "see-pee two-thirty-three" |
| E119Q | "E one-nineteen Q" |
| R146 / R146A | "R one-forty-six" / "R one-forty-six A" |
| Y129 | "Y one-twenty-nine" |
| Å / angstrom | "angstrom" |
| pN | "piconewton" |
| ClpX | "clip-ex" |
| MD | "molecular dynamics" |
| MSA | "em-ess-ay" |
| Boltz-2 / OpenFold3 | "bolts two" / "open-fold three" |

## Delivery reminders
- Slow on numbers; the visuals hold while you say them.
- The two beats that win trust are **v08** (the honest null) and **v14** (cross-checked, not validated). Do not throw them away.
- Never say "pTM" as evidence — the whole point is we *don't* steer by it.
