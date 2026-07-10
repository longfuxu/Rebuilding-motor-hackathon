#!/usr/bin/env python3
"""High-quality FREE voiceover via Microsoft Edge neural TTS (edge-tts).
No API key, no credits — genuinely neural voices, near-commercial quality.
Writes audio/v01.mp3 … v14.mp3.

  pip install edge-tts
  python generate_voiceover_edge.py                    # default voice
  VOICE="en-US-BrianMultilingualNeural" python generate_voiceover_edge.py

Good documentary voices:
  en-US-AndrewMultilingualNeural  (warm, natural male — default)
  en-US-BrianMultilingualNeural   (deep, calm male)
  en-US-GuyNeural                 (news/narrator male)
  en-US-AriaNeural / en-US-JennyNeural (female)
  en-GB-RyanNeural                (British male)
"""
import os, asyncio, pathlib
import edge_tts
from narration_lines import LINES

KIT = pathlib.Path(__file__).resolve().parent
OUT = KIT / "audio"; OUT.mkdir(exist_ok=True)

VOICE = os.environ.get("VOICE", "en-US-AndrewMultilingualNeural")
RATE  = os.environ.get("RATE", "-4%")     # slightly slower = more authoritative
PITCH = os.environ.get("PITCH", "+0Hz")

async def one(vid, text):
    tts = edge_tts.Communicate(text, VOICE, rate=RATE, pitch=PITCH)
    await tts.save(str(OUT / f"{vid}.mp3"))
    kb = (OUT / f"{vid}.mp3").stat().st_size // 1024
    print(f"[{vid}] {kb} KB · {VOICE}")

async def main():
    ids = [a for a in os.sys.argv[1:] if a in LINES] or list(LINES)
    for vid in ids:
        await one(vid, LINES[vid])
    print(f"\n{len(ids)} MP3s written to {OUT}  (voice: {VOICE})")

if __name__ == "__main__":
    asyncio.run(main())
