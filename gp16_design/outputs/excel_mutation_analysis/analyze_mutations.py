#!/usr/bin/env python3.13
"""
phi29 gp16 bulk-biochemistry mutation dataset  ->  cross-check against our
single-chain design proxy metrics (M2/M3/M4/M5 functional residues).

Answers 3 questions:
  Q1  Do biochem-important residues OVERLAP our proxy-flagged residues?
  Q2  Does CP233 avoid biochem-critical residues while cp285/cp297 do not?
  Q3  Tolerant vs lethal residue lists = prior for directed evolution.

Input : phi29_mapped_residues.xlsx (Sheet1)
Output: cleaned_residue_phenotype.csv  +  numbers printed for REPORT.md
"""
import re
import pandas as pd
from collections import defaultdict

XLSX = ('/Users/longfu/Library/CloudStorage/Dropbox/2-Postdoc Research/'
        '02_Phage Protocol Files/phi29_mapped_residues.xlsx')
OUT  = ('/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/'
        'export-gp16-results/gp16_design/outputs/excel_mutation_analysis/'
        'cleaned_residue_phenotype.csv')

# ---------------------------------------------------------------- proxy sets
def expand(*rngs):
    s = set()
    for r in rngs:
        if isinstance(r, tuple):
            s.update(range(r[0], r[1] + 1))
        else:
            s.add(r)
    return s

DNA_CONTACT = {55,56,57,59,60,82,83,98,99,100,125,126,127,128,129,130,292,293,297,330}
PRNA_CONTACT= {10,11,13,14,15,16,17,37,40,41,44,45,148,149,150,152,182,
               237,238,239,240,254,267,268,269,270,271,272,273}
WALKER_A    = expand((24,31))          # P-loop
WALKER_B    = {118,119}                # D118 / E119
Y129        = {129}                    # force-transmission pivot
R146        = {146}                    # trans arginine finger (M2 soul)
CP_SITES    = {233,285,297}            # circular-permutation cut points

def proxy_tags(n):
    t = []
    if n in DNA_CONTACT:  t.append('DNA-contact')
    if n in PRNA_CONTACT: t.append('pRNA-contact')
    if n in WALKER_A:     t.append('Walker-A')
    if n in WALKER_B:     t.append('Walker-B')
    if n in Y129:         t.append('Y129-force')
    if n in R146:         t.append('R146-transfinger')
    if n in CP_SITES:     t.append('CP-site')
    return t

PROXY_ALL = DNA_CONTACT | PRNA_CONTACT | WALKER_A | WALKER_B | Y129 | R146 | CP_SITES

# ---------------------------------------------------------------- parsing
AA3 = {'Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu',
       'Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val'}
AA1 = set('ACDEFGHIKLMNPQRSTVWY')

def res_num_from_label(x):
    """'Tyr297' -> (297,'Tyr')"""
    if not isinstance(x, str): return None
    m = re.match(r'([A-Z][a-z]{2})(\d+)', x.strip())
    if m and m.group(1) in AA3:
        return int(m.group(2)), m.group(1)
    return None

def muts_from_made(x):
    """'K233A-R234A' -> [(233,'K','A'),(234,'R','A')]; 'TO DO: F6A' -> [(6,'F','A')]"""
    if not isinstance(x, str): return []
    out = []
    for m in re.finditer(r'([A-Z])(\d+)([A-Z])', x):
        wt, num, mt = m.group(1), int(m.group(2)), m.group(3)
        if wt in AA1 and mt in AA1:
            out.append((num, wt, mt))
    # single form 'F6A' already caught; also 'F6' with no target (rare) ignore
    return out

def to_frac(txt):
    """normalize a packaging/ATPase phenotype string -> fraction of Wt (0..1) or None"""
    if txt is None: return None
    if isinstance(txt,float) and pd.isna(txt): return None      # <-- NaN != lethal
    if isinstance(txt,(int,float)):
        # ClaI column uses 0 / 1 / 0.2 style
        return float(txt)
    s = str(txt).lower().strip()
    if s in ('nan','','??'): return None
    if 'insoluble' in s: return None                            # folding, not a fraction
    if 'trace' in s or s.startswith('0%'): return 0.0
    # ranges like '60-70%'
    r = re.search(r'(\d+)\s*-\s*(\d+)\s*%', s)
    if r: return (int(r.group(1))+int(r.group(2)))/200.0
    p = re.search(r'(\d+(?:\.\d+)?)\s*%', s)
    if p: return float(p.group(1))/100.0
    if s in ('0','0.0'): return 0.0
    return None

# ---------------------------------------------------------------- load
df = pd.read_excel(XLSX, sheet_name='Sheet1')

# residue label carried on some rows
residue_of_row = []
last = None
for v in df['residue']:
    r = res_num_from_label(v)
    if r: last = r
    residue_of_row.append(last)

# accumulate per-residue info
info = defaultdict(lambda: {
    'aa':'', 'interaction':set(), 'predicted':set(),
    'measured':[],      # clean single-mutant readouts
    'epistatic':[],     # double-mutant / background readouts (not single-residue tolerance)
})

cur_interaction = None
for i, row in df.iterrows():
    if isinstance(row['interaction'], str) and row['interaction'].strip():
        cur_interaction = row['interaction'].strip()
    # predicted phenotype attaches to the residue label on this row
    lab = res_num_from_label(row['residue'])
    if lab:
        n, aa = lab
        info[n]['aa'] = aa
        if cur_interaction: info[n]['interaction'].add(cur_interaction)
        pp = row['predicted mutant phenotype']
        if isinstance(pp,str) and pp.strip():
            info[n]['predicted'].add(pp.strip())
    # measured mutants on this row
    made_txt = str(row['mutant made'])
    muts = muts_from_made(row['mutant made'])
    # a cell with >1 mutation, or an explicit "background", is an epistatic combo
    is_bg = ('background' in made_txt.lower()) or (len(muts) > 1)
    for (n, wt, mt) in muts:
        info[n]['aa'] = info[n]['aa'] or wt
        if cur_interaction: info[n]['interaction'].add(cur_interaction)
        raw_pack = row['bulk packaging phenotype']
        raw_clai = row['ClaI fragment packaging']
        raw_atp  = row['bulk prohead/gp16 ATPase']
        pack = to_frac(raw_pack)
        clai = to_frac(raw_clai)
        atp  = to_frac(raw_atp)
        # combined packaging readout = prefer bulk, else ClaI
        pk = pack if pack is not None else clai
        entry = dict(mut=f'{wt}{n}{mt}', pack=pk, atp=atp,
                     raw_pack=raw_pack, raw_clai=raw_clai, raw_atp=raw_atp,
                     context=made_txt.strip())
        if is_bg:
            info[n]['epistatic'].append(entry)   # e.g. Q222x in K105A background
        else:
            info[n]['measured'].append(entry)

# ---------------------------------------------------------------- classify
# conservative-substitution groups (side-chain chemistry preserved)
CONS_GROUPS = [set('KR'), set('DE'), set('ST'), set('FYW'),
               set('ILVM'), set('NQ'), set('AG')]
def is_conservative(wt, mt):
    return any(wt in g and mt in g for g in CONS_GROUPS)

def classify(rec):
    """tolerance class from measured mutants; packaging preferred over ATPase."""
    ms = [m for m in rec['measured'] if m['pack'] is not None or m['atp'] is not None]
    if not ms:
        made = bool(rec['measured'])
        if rec['predicted']: return 'PREDICTED-ONLY'
        return 'MADE-NO-DATA' if made else 'UNTESTED'
    # per-mutant functional readout: packaging if present else ATPase
    vals = []
    for m in ms:
        v = m['pack'] if m['pack'] is not None else m['atp']
        if v is not None:
            vals.append((v, m['mut']))
    if not vals: return 'MADE-NO-DATA'
    best = max(v for v, _ in vals)
    if   best >= 0.70: return 'TOLERANT'
    elif best >= 0.30: return 'PARTIAL'
    elif best >= 0.05: return 'SEVERE'
    else:              return 'LETHAL'

def conservative_note(rec):
    """flag residues tolerated ONLY under a conservative substitution
    (chemistry-critical: a non-conservative sub breaks it, a conservative one rescues)."""
    data = []
    for m in rec['measured']:
        v = m['pack'] if m['pack'] is not None else m['atp']
        if v is None: continue
        data.append((v, m['mut'][0], m['mut'][-1]))
    ala = next((v for v, wt, mt in data if mt == 'A'), None)
    if not data:
        return '', ala
    best = max(v for v, _, _ in data)
    best_is_cons = any(v >= 0.70 and is_conservative(wt, mt) for v, wt, mt in data)
    noncons_break = any(v < 0.30 and not is_conservative(wt, mt) for v, wt, mt in data)
    note = []
    if best >= 0.70 and best_is_cons and noncons_break:
        note.append('chemistry-critical(only conservative sub tolerated)')
    return ';'.join(note), ala

# axis = which functional machine the residue belongs to
# force / portal-CTD grip + pRNA-CTD dynamic contacts (bind DNA/pRNA when CTD rotates)
FORCE_PORTAL = {175,176,177,178, 233,234,236, 294,300,301,327,328}
GRIP_SET = DNA_CONTACT | PRNA_CONTACT | FORCE_PORTAL | Y129
def axis_of(n):
    return 'grip' if n in GRIP_SET else 'atp-active-site'

rows = []
for n in sorted(info):
    rec = info[n]
    cls = classify(rec)
    tags = proxy_tags(n)
    note, ala = conservative_note(rec)
    packs = [m['pack'] for m in rec['measured'] if m['pack'] is not None]
    atps  = [m['atp']  for m in rec['measured'] if m['atp'] is not None]
    insoluble = any('insoluble' in str(m['raw_pack']).lower() or
                    'insoluble' in m['mut'] for m in rec['measured'])
    rows.append(dict(
        residue=n,
        aa=rec['aa'],
        interaction=';'.join(sorted(rec['interaction'])),
        axis=axis_of(n),
        proxy_flagged=bool(tags),
        proxy_tags=';'.join(tags),
        n_mutants=len(rec['measured']),
        mutants=';'.join(m['mut'] for m in rec['measured']),
        best_pack=round(max(packs),3) if packs else '',
        worst_pack=round(min(packs),3) if packs else '',
        ala_frac=round(ala,3) if ala is not None else '',
        best_atp=round(max(atps),3) if atps else '',
        worst_atp=round(min(atps),3) if atps else '',
        tolerance_class=cls,
        chem_note=note,
        insoluble=insoluble,
        epistatic=';'.join(f"{m['mut']}[{m['context']}]={m['pack']}"
                           for m in rec['epistatic']),
        predicted='; '.join(sorted(rec['predicted'])),
    ))

out = pd.DataFrame(rows).sort_values('residue')
out.to_csv(OUT, index=False)
print('WROTE', OUT, out.shape)

# ---------------------------------------------------------------- Q1 overlap
print('\n=========== Q1  OVERLAP (proxy vs biochem) ===========')
# measured residues with a clear defect/tolerance call
measured = out[out['tolerance_class'].isin(
    ['LETHAL','SEVERE','PARTIAL','TOLERANT'])].copy()
defects  = measured[measured['tolerance_class'].isin(['LETHAL','SEVERE','PARTIAL'])]
print(f'residues measured with a call : {len(measured)}')
print(f'  of which DEFECT (L/S/P)     : {len(defects)}')
print(f'  of which TOLERANT           : {(measured.tolerance_class=="TOLERANT").sum()}')

flagged_def = defects[defects.proxy_flagged]
print(f'\nDEFECT residues proxy-flagged (exact): {len(flagged_def)}/{len(defects)}'
      f' = {100*len(flagged_def)/len(defects):.0f}%')
print('  flagged:', list(flagged_def.residue))
print('  MISSED :', list(defects[~defects.proxy_flagged].residue))

# within +-1 residue of any proxy residue
def near(n, d=1):
    return any((n+k) in PROXY_ALL for k in range(-d, d+1))
defects = defects.copy()
defects['near_proxy'] = defects.residue.apply(near)
print(f'\nDEFECT residues within +-1 of a proxy residue:'
      f' {defects.near_proxy.sum()}/{len(defects)}'
      f' = {100*defects.near_proxy.sum()/len(defects):.0f}%')
print('  still missed (>1 away):', list(defects[~defects.near_proxy].residue))

# split by functional axis (assigned by residue membership, see axis_of)
for ax in ('grip','atp-active-site'):
    sub = defects[defects.axis==ax]
    fl = sub.proxy_flagged.sum()
    nr = sub.near_proxy.sum()
    print(f'\n  axis={ax}: n={len(sub)}  exact-flag={fl} ({100*fl/max(len(sub),1):.0f}%)'
          f'  within1={nr} ({100*nr/max(len(sub),1):.0f}%)')
    print('    residues:', list(sub.residue))

# ---------------------------------------------------------------- Q2 CP sites
print('\n=========== Q2  CP SITE vs CRITICAL RESIDUES ===========')
crit = set(defects[defects.tolerance_class.isin(['LETHAL','SEVERE'])].residue)
print('measured LETHAL/SEVERE residues:', sorted(crit))
for cp in sorted(CP_SITES):
    rec = info.get(cp)
    cls = out[out.residue==cp]['tolerance_class'].values
    cls = cls[0] if len(cls) else 'UNTESTED'
    d_lo = min((abs(cp-c) for c in crit), default=None)
    nearest = min(crit, key=lambda c: abs(cp-c)) if crit else None
    incontact = cp in DNA_CONTACT or cp in PRNA_CONTACT
    print(f'  CP{cp}: cut-residue class={cls:14s} '
          f'in_contact_face={incontact}  nearest_critical={nearest} (|d|={d_lo})')

# ---------------------------------------------------------------- Q3 lists
print('\n=========== Q3  DIRECTED-EVOLUTION PRIOR ===========')
tol = out[out.tolerance_class=='TOLERANT']
leth= out[out.tolerance_class=='LETHAL']
sev = out[out.tolerance_class=='SEVERE']
par = out[out.tolerance_class=='PARTIAL']
print('TOLERANT (safe to mutate):',
      [f"{r.aa}{r.residue}" for r in tol.itertuples()])
print('LETHAL   (do not touch)  :',
      [f"{r.aa}{r.residue}" for r in leth.itertuples()])
print('SEVERE   (near-lethal)   :',
      [f"{r.aa}{r.residue}" for r in sev.itertuples()])
print('PARTIAL  (tunable)       :',
      [f"{r.aa}{r.residue}" for r in par.itertuples()])

# chemistry-constrained residues (only conservative sub tolerated)
chem = out[out.chem_note.astype(str).str.contains('chemistry-critical')]
print('\nCHEMISTRY-CONSTRAINED (conservative sub only):',
      [f"{r.aa}{r.residue}" for r in chem.itertuples()])

# proxy DNA-contact residues split into load-bearing vs redundant
print('\n--- DNA-contact residues that were TESTED (load-bearing vs redundant) ---')
for r in out.itertuples():
    if r.residue in DNA_CONTACT and r.tolerance_class in (
        'LETHAL','SEVERE','PARTIAL','TOLERANT'):
        print(f'  {r.aa}{r.residue}: {r.tolerance_class:9s} ({r.mutants})')

# ---------------------------------------------------------------- Q1b predicted
print('\n=========== Q1b  PREDICTED-important residues vs proxy ===========')
pred = out[out.predicted.astype(str).str.len() > 0]
pin  = pred[pred.proxy_flagged]
print(f'residues with a PREDICTED phenotype: {len(pred)}')
print(f'  proxy-flagged: {len(pin)}/{len(pred)} = {100*len(pin)/len(pred):.0f}%')
print('  flagged   :', [f"{r.aa}{r.residue}" for r in pin.itertuples()])
print('  NOT flagged:', [f"{r.aa}{r.residue}" for r in pred[~pred.proxy_flagged].itertuples()])

# ---------------------------------------------------------------- proxy gaps
print('\n=========== PROXY GAPS (biochem-critical, our proxy misses) ===========')
crit_all = out[out.tolerance_class.isin(['LETHAL','SEVERE'])]
gap = crit_all[~crit_all.proxy_flagged]
for r in gap.itertuples():
    print(f'  {r.aa}{r.residue}  axis={r.axis}  interaction={r.interaction}')
print('\nNOTE: JSON dna_contact_res omits 58 though METRICS.md writes "55-60";'
      ' measured E58A is SEVERE (glu-switch). Range vs enumerated-set mismatch.')
