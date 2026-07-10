#!/usr/bin/env python3
"""Voiceover via ElevenLabs — highest-quality model, one MP3 per slide.

The API key is read at runtime from the environment or a .env file; it is never
printed or written to any output. Lookup order for ELEVENLABS_API_KEY:
  1. process environment   2. <repo>/.env   3. ~/.hermes/.env   4. ~/.env

  python generate_voiceover.py                 # best-quality defaults, all lines
  VOICE_ID=<id> python generate_voiceover.py   # pick a specific ElevenLabs voice
  MODEL_ID=eleven_v3 python generate_voiceover.py   # most expressive (if enabled)

NOTE: needs ElevenLabs credits. A near-commercial FREE alternative that needs no
key is generate_voiceover_edge.py (Microsoft neural TTS).
"""
import os, sys, pathlib, re
from narration_lines import LINES

KIT = pathlib.Path(__file__).resolve().parent
OUT = KIT / "audio"; OUT.mkdir(exist_ok=True)

def load_key():
    # case-insensitive env lookup
    for name, val in os.environ.items():
        if name.upper() == "ELEVENLABS_API_KEY" and val.strip():
            return val.strip()
    for f in [pathlib.Path("/Users/longfu/Developer/claude-science-hackthon/.env"),
              pathlib.Path.home() / ".hermes" / ".env",
              pathlib.Path.home() / ".env"]:
        if f.is_file():
            for line in f.read_text().splitlines():
                m = re.match(r'\s*(?:export\s+)?ELEVENLABS_API_KEY\s*=\s*(.+)\s*$', line, re.IGNORECASE)
                if m:
                    return m.group(1).strip().strip('"').strip("'")
    return None

# highest-quality defaults
VOICE_ID = os.environ.get("VOICE_ID", "pNInz6obpgDQGcFmaJgB")     # Adam (calm narrator)
MODEL_ID = os.environ.get("MODEL_ID", "eleven_multilingual_v2")   # best quality (not turbo/flash)
OUT_FMT  = os.environ.get("OUT_FMT",  "mp3_44100_192")            # highest MP3 bitrate (needs paid tier)
VOICE_SETTINGS = {"stability": 0.45, "similarity_boost": 0.80, "style": 0.15, "use_speaker_boost": True}

def main():
    key = load_key()
    if not key:
        sys.exit("ELEVENLABS_API_KEY not found. Add  ELEVENLABS_API_KEY=...  to .env and re-run, "
                 "or use the free generate_voiceover_edge.py instead.")
    import requests
    url = f"https://api.elevenlabs.io/v1/text-to-speech/{VOICE_ID}?output_format={OUT_FMT}"
    head = {"xi-api-key": key, "Content-Type": "application/json", "Accept": "audio/mpeg"}
    ids = [a for a in sys.argv[1:] if a in LINES] or list(LINES)
    ok = 0
    for vid in ids:
        r = requests.post(url, headers=head,
                          json={"text": LINES[vid], "model_id": MODEL_ID, "voice_settings": VOICE_SETTINGS},
                          timeout=180)
        if r.status_code != 200:
            # if 192kbps needs a paid tier, retry once at 128
            r2 = requests.post(f"https://api.elevenlabs.io/v1/text-to-speech/{VOICE_ID}?output_format=mp3_44100_128",
                               headers=head, json={"text": LINES[vid], "model_id": MODEL_ID, "voice_settings": VOICE_SETTINGS}, timeout=180)
            if r2.status_code == 200:
                r = r2
            else:
                print(f"[{vid}] HTTP {r.status_code}: {r.text[:180]}")
                continue
        (OUT / f"{vid}.mp3").write_bytes(r.content)
        ok += 1; print(f"[{vid}] wrote {len(r.content)//1024} KB")
    print(f"\n{ok}/{len(ids)} MP3s written to {OUT}")

if __name__ == "__main__":
    main()
