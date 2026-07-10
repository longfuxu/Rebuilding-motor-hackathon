#!/usr/bin/env python3
"""Render the demo deck to 1920x1080 slide PNGs (Pillow) and assemble a narrated
MP4 with ffmpeg, each slide held for the length of its voiceover clip.

  python build_video.py            # full 14-slide video + 8-slide spine cut
Outputs: video/gp16_demo_full.mp4, video/gp16_demo_spine.mp4, video_slides/*.png
Audio: uses audio/vNN.mp3 (run generate_voiceover*.py first).
"""
import os, subprocess, pathlib
from PIL import Image, ImageDraw, ImageFont

KIT = pathlib.Path(__file__).resolve().parent
FIG = KIT / "figs"; AUD = KIT / "audio"
SL  = KIT / "video_slides"; SL.mkdir(exist_ok=True)
VID = KIT / "video"; VID.mkdir(exist_ok=True)
W, H = 1920, 1080

# palette (RGB)
GROUND=(11,15,20); PANEL=(18,26,32); INK=(238,244,243); DIM=(159,178,176); FAINT=(103,128,125)
TEAL=(55,211,197); AMBER=(230,171,58); CORAL=(236,106,134); WHITE=(246,248,247); LINE=(36,51,48)

DISP="/System/Library/Fonts/Supplemental/Arial.ttf"
DISPB="/System/Library/Fonts/Supplemental/Arial Bold.ttf"
MONO="/System/Library/Fonts/Menlo.ttc"
_fc={}
def font(size, bold=False, mono=False):
    key=(size,bold,mono)
    if key not in _fc:
        path = MONO if mono else (DISPB if bold else DISP)
        _fc[key]=ImageFont.truetype(path, size)
    return _fc[key]

def tw(d, s, f): return d.textlength(s, font=f)

def draw_runs(d, x, y, runs, f):
    for s,c in runs:
        d.text((x,y), s, font=f, fill=c); x += tw(d,s,f)
    return x

def wrap_runs(d, runs, f, maxw):
    toks=[]
    for s,c in runs:
        for w in s.split(' '):
            if w!='': toks.append((w,c))
    lines=[]; cur=[]; curw=0; sp=tw(d,' ',f)
    for w,c in toks:
        ww=tw(d,w,f); add=ww+(sp if cur else 0)
        if cur and curw+add>maxw:
            lines.append(cur); cur=[(w,c)]; curw=ww
        else:
            cur.append((w,c)); curw+=add
    if cur: lines.append(cur)
    return lines

def draw_wrapped(d, x0, y, runs, f, maxw, lh=None):
    lh = lh or int(f.size*1.16)
    for line in wrap_runs(d, runs, f, maxw):
        x=x0; sp=tw(d,' ',f)
        for i,(w,c) in enumerate(line):
            if i: x+=sp
            d.text((x,y), w, font=f, fill=c); x+=tw(d,w,f)
        y+=lh
    return y

def kicker(d, text):
    f=font(26, mono=True)
    draw_runs(d, 96, 66, [("—  ",TEAL),(text.upper(),TEAL)], f)

def footer(d, right=""):
    d.line([(96,986),(1824,986)], fill=LINE, width=2)
    f=font(20, mono=True)
    draw_runs(d, 96, 1002, [("Built with ",FAINT),("Claude Science",TEAL),(" + ",FAINT),("Claude Code",TEAL)], f)
    if right:
        d.text((1824-tw(d,right,f),1002), right, font=f, fill=FAINT)

def chips(d, items, y0=858):
    f=font(25, mono=True); x0=96; x=x0; y=y0; gap=42; r=6; lh=int(f.size*1.75); maxw=1728
    for i,c in enumerate(items):
        cw=r*2+12+tw(d,c,f)
        if x>x0 and x+cw>x0+maxw: x=x0; y+=lh
        cy=y+f.size//2
        d.ellipse([x,cy-r,x+2*r,cy+r], fill=[TEAL,AMBER,CORAL][i%3])
        d.text((x+2*r+10,y), c, font=f, fill=INK)
        x+=cw+gap

def figure_card(im, d, figname, cx=96, cy=300, cw=1728, ch=536, pad=26, crop=None):
    d.rounded_rectangle([cx,cy,cx+cw,cy+ch], radius=18, fill=WHITE)
    fig=Image.open(FIG/figname).convert("RGB")
    if crop:  # (l,t,r,b) as fractions
        l,t,r,b=crop; fw,fh=fig.size
        fig=fig.crop((int(l*fw),int(t*fh),int(r*fw),int(b*fh)))
    mw,mh=cw-2*pad, ch-2*pad; a=fig.width/fig.height
    if mw/mh>a: hh=mh; ww=int(mh*a)
    else: ww=mw; hh=int(mw/a)
    fig=fig.resize((ww,hh), Image.LANCZOS)
    im.paste(fig, (cx+(cw-ww)//2, cy+(ch-hh)//2))

def base():
    im=Image.new("RGB",(W,H),GROUND); return im, ImageDraw.Draw(im)

# ---------------- slide renderers ----------------
def render(i, s):
    im,d=base(); k=s["kind"]
    if k=="title":
        kicker(d, s["kicker"])
        y=draw_wrapped(d, 96, 320, s["title"], font(96,bold=True), 1660, lh=110)
        y=draw_wrapped(d, 96, y+30, [(s["lead"],DIM)], font(30), 1360, lh=44)
        draw_runs(d, 96, 880, s["meta"], font(24))
        footer(d, s.get("foot",""))
    elif k=="close":
        kicker(d, s["kicker"])
        y=draw_wrapped(d, 96, 300, s["title"], font(80,bold=True), 1720, lh=94)
        y=draw_wrapped(d, 96, y+34, [(s["lead"],DIM)], font(29), 1500, lh=42)
        draw_runs(d, 96, 880, s["meta"], font(24))
        footer(d, s.get("foot",""))
    elif k=="text":
        kicker(d, s["kicker"])
        y=draw_wrapped(d, 96, 150, s["head"], font(60,bold=True), 1720, lh=72)
        y=draw_wrapped(d, 96, y+40, [(s["body"],INK)], font(34), 1500, lh=50)
        if s.get("accent"):
            draw_wrapped(d, 96, y+30, [s["accent"]], font(28,mono=True), 1600, lh=40)
        if s.get("chips"): chips(d, s["chips"])
        footer(d, s.get("foot",""))
    elif k=="bullets":
        kicker(d, s["kicker"])
        draw_wrapped(d, 96, 150, s["head"], font(58,bold=True), 1720, lh=70)
        y=290
        for hd,bd in s["bullets"]:
            ty=y+9; d.polygon([(96,ty),(96,ty+24),(116,ty+12)], fill=TEAL)  # marker
            d.text((134,y), hd, font=font(30,bold=True), fill=INK)
            y=draw_wrapped(d, 134, y+46, [(bd,DIM)], font(27), 1620, lh=38)+18
        if s.get("chips"): chips(d, s["chips"])
        footer(d, s.get("foot",""))
    else:  # figure
        kicker(d, s["kicker"])
        draw_wrapped(d, 96, 138, s["head"], font(58,bold=True), 1720, lh=70)
        figure_card(im, d, s["fig"], crop=s.get("crop"))
        chips(d, s["chips"])
        footer(d, s.get("foot",""))
    out=SL/f"slide_{i:02d}.png"; im.save(out); return out

SLIDES=[
 {"kind":"title","kicker":"Built with Claude · Life Sciences — Research track",
  "title":[("Rebuilding a Molecular Motor ",INK),("to Understand It",TEAL)],
  "lead":"A single-chain, genetically-addressable φ29 gp16 DNA-packaging ring — designed and cross-validated end-to-end.",
  "meta":[("Longfu Xu",INK),("   ·   Bustamante Lab, UC Berkeley / HHMI",DIM)],
  "foot":"gp16 · PDB 7JQQ · UniProt P11014"},
 {"kind":"text","kicker":"The problem",
  "head":[("Nature's strongest motors are rings of ",INK),("identical",AMBER),(" parts.",INK)],
  "body":"φ29 gp16 packs DNA at ~57 pN as a homo-pentamer. Its deepest questions are per-subunit: which seat is the special, regulatory subunit? Is ATP hydrolysis sequential or concerted? In a ring of identical parts you cannot put a defect at a known seat — every measurement averages over an unknown mixture.",
  "chips":["5 identical subunits","4 workers + 1 special (Chistol 2012)","status quo: stochastic doping"],
  "foot":"the special-subunit problem"},
 {"kind":"text","kicker":"The move",
  "head":[("Rebuild the ring as ",INK),("one addressable chain.",TEAL)],
  "body":"Fuse five protomers into a single gene: a mutation encoded in copy k lands at ring seat k — deterministically. The trick that made single-chain ClpX interpretable (Martin, Baker & Sauer 2005). Never reported for gp16.",
  "accent":("Understanding by rebuilding — build it better than evolution did, and you've proven you understand it.",AMBER),
  "foot":"single-chain, position-addressable"},
 {"kind":"figure","kicker":"Result · the lead is real",
  "head":[("cp233",TEAL),(" — one 1,750-aa chain that folds into the native ring.",INK)],
  "fig":"02_cp233_vs_native.png",
  "chips":["3 predictors agree 5/5 (Boltz-2 · OpenFold3 · AlphaFold3)","subunit RMSD 1.80 Å · TM 0.94","all 5 copies fall on the native ring"],
  "foot":"tiled-MSA fold · steered by M1/M2, never pTM"},
 {"kind":"figure","kicker":"Result · physics, not just predictors",
  "head":[("MD keeps the design closed — while the ",INK),("ATP-bound native opens.",AMBER)],
  "fig":"04_md_per_interface.png",
  "chips":["apo ring stays closed","7JQQ opens 2 seams","design stays 5/5 engaged","interface −192 vs −174 kcal/mol"],
  "foot":"OpenMM GBSA · predictor-independent"},
 {"kind":"figure","kicker":"Result · function is a stricter filter",
  "head":[("It threads DNA. Look-alikes that merely ",INK),("“close” don’t.",CORAL)],
  "fig":"03_cp233_threads_dna.png",
  "chips":["cp233 channel 20.8 Å ≈ native 20.6 Å","cp285/cp297 close but can't thread","closes ≠ works"],
  "foot":"M3 · DNA-channel competence"},
 {"kind":"figure","kicker":"The payoff · genetic addressability",
  "head":[("A catalytically-dead seat at a ",INK),("chosen position.",TEAL)],
  "fig":"01_deadseat_addressable.png","crop":(0,0,1,0.54),
  "chips":["E119Q at 1 or all 5 seats → ring stays 5/5 closed","structurally-silent: Mg not displaced","R146A localizes the defect"],
  "foot":"apo + ATP·Mg · 3-predictor cross-checked"},
 {"kind":"figure","kicker":"Result · and an honest null",
  "head":[("On gp16 the design out-couples native — but it ",INK),("doesn't generalize.",CORAL)],
  "fig":"05_clpx_vs_gp16_coupling.png",
  "chips":["gp16 coupling 1.79×","Y129 hub 4.6× more central","ClpX ≈ parity — we report the null"],
  "foot":"advantage is topological, not covalency"},
 {"kind":"figure","kicker":"The framework · it generalizes",
  "head":[("13 ring motors collapse to ",INK),("one predictive build rule.",TEAL)],
  "fig":"06_buildability_atlas.png",
  "chips":["terminus fouls channel → circular permutation","otherwise → direct fusion","3/3 on validated anchors: gp16 · ClpX · gp17"],
  "foot":"from structure alone, predict the build"},
 {"kind":"figure","kicker":"Credibility · validated on a known answer",
  "head":[("Run blind on ClpX, the framework catches the ",INK),("dead couplers.",CORAL)],
  "fig":"11_clpx_retro_validation.png",
  "chips":["good WT · 6/6 interfaces engaged","R307A / R307E · 0/6 (caught)","M2 separates functional from broken"],
  "foot":"retrospective external validation"},
 {"kind":"figure","kicker":"Method · why it was hard",
  "head":[("A function-first metric hierarchy — ",INK),("never global pTM.",CORAL)],
  "fig":"10_msa_ladder.png",
  "chips":["M1 ring · M2 R146 · M3 DNA · M4 pRNA · M5 ATP","MSA-tiling ≈ 2× pTM made the 1,750-aa ring foldable"],
  "foot":"steer by geometry, not confidence"},
 {"kind":"figure","kicker":"What it unlocks · follow-on program → follow-on",
  "head":[("Set the ATP pattern, read the force — ",INK),("sequential vs concerted.",AMBER)],
  "fig":"12_mechanism_atp_occupancy.png",
  "chips":["opening is progressive & cooperative","which subunits fire sets the state","addressable ring + tweezers = the experiment"],
  "foot":"the design → measurement bridge"},
 {"kind":"bullets","kicker":"Built with Claude",
  "head":[("Claude as a ",INK),("co-scientist,",TEAL),(" not a wrapper.",INK)],
  "bullets":[
   ("Claude Code = orchestrator.","Drove every GPU job (RFdiffusion, MD, NIM folding) through reliable hard-timeout harnesses; fanned out a dozen parallel agents."),
   ("A 19-sub-agent workflow.","Score → adversarially verify → synthesize every AlphaFold3 fold."),
   ("Adversarial verification changed the science.","It caught handedness false-negatives in M2, and surfaced the ClpX null and the pLDDT-trap honestly."),
   ("Claude Science","ran the systematic folding & validation campaigns and rendered the structures.")],
  "chips":["< $10 total compute","128 commits · fully scripted","packaged as a reusable protein-design skill"],
  "foot":"actor–critic, at scale"},
 {"kind":"close","kicker":"The closed loop",
  "title":[("Predict the build. Rank the coupling. ",INK),("Design stronger. Then measure.",TEAL)],
  "lead":"Cross-checked, not yet wet-lab-validated — every prediction is a falsifiable hypothesis for single-molecule assays. That honesty is the point.",
  "meta":[("Longfu Xu",INK),("   ·   a designed motor we can both predict and measure.",DIM)],
  "foot":"thank you"},
]

def dur(mp3):
    r=subprocess.run(["ffprobe","-v","error","-show_entries","format=duration",
                      "-of","default=nw=1:nk=1",str(mp3)],capture_output=True,text=True)
    return float(r.stdout.strip())

def clip(i):
    png=SL/f"slide_{i:02d}.png"; mp3=AUD/f"v{i:02d}.mp3"; out=VID/f"clip_{i:02d}.mp4"
    subprocess.run(["ffmpeg","-y","-loop","1","-i",str(png),"-i",str(mp3),
        "-c:v","libx264","-preset","medium","-tune","stillimage","-pix_fmt","yuv420p",
        "-r","30","-c:a","aac","-b:a","192k","-shortest",str(out)],
        check=True, capture_output=True)
    return out

def concat(clips, out):
    lst=VID/"_list.txt"
    lst.write_text("".join(f"file '{c}'\n" for c in clips))
    subprocess.run(["ffmpeg","-y","-f","concat","-safe","0","-i",str(lst),"-c","copy",str(out)],
                   check=True, capture_output=True)
    lst.unlink()

if __name__=="__main__":
    print("rendering slides…")
    for i,s in enumerate(SLIDES,1):
        render(i,s); print(f"  slide {i:02d}  {s['kind']}")
    print("building clips…")
    clips=[clip(i) for i in range(1,15)]
    concat(clips, VID/"gp16_demo_full.mp4")
    spine=[VID/f"clip_{i:02d}.mp4" for i in (1,2,3,4,7,9,13,14)]
    concat(spine, VID/"gp16_demo_spine.mp4")
    for c in clips: c.unlink()
    for f in ["gp16_demo_full.mp4","gp16_demo_spine.mp4"]:
        d=dur(VID/f); print(f"{f}: {int(d//60)}:{int(d%60):02d}")
    print("done ->", VID)
