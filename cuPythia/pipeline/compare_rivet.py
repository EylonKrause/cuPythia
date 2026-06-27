#!/usr/bin/env python3
# Compare a cuPythia Rivet histogram (Histo1D) to the ALEPH reference (Scatter2D) for one observable.
# Usage: compare_rivet.py <mc.yoda> <ref.yoda> <histo-path-suffix>   e.g. d54-x01-y01 (thrust @91.2)
import sys, math

def grab(path, beginpat, kind):
    out, f = [], False
    for ln in open(path, encoding='utf-8', errors='ignore'):
        if ln.startswith('BEGIN YODA') and beginpat in ln: f=True; continue
        if f and ln.startswith('END YODA'): break
        if f and ln and ln[0].isdigit(): out.append(ln.split())
        if f and kind=='histo' and (ln.startswith('Total') or ln.startswith('Underflow') or ln.startswith('Overflow')):
            continue
    return out

mc_f, ref_f, suf = sys.argv[1], sys.argv[2], sys.argv[3]
# MC Histo1D bins: xlow xhigh sumw sumw2 sumwx sumwx2 nEntries -> density = sumw/width
mc=[]
for r in grab(mc_f, '/ALEPH_2004_S5765862/'+suf, 'histo'):
    if len(r)>=3:
        xlo,xhi,sumw=float(r[0]),float(r[1]),float(r[2])
        if xhi>xlo: mc.append((0.5*(xlo+xhi), xlo, xhi, sumw/(xhi-xlo)))
# REF Scatter2D points: xval xerr- xerr+ yval yerr- yerr+
ref=[]
for r in grab(ref_f, '/REF/ALEPH_2004_S5765862/'+suf, 'scatter'):
    if len(r)>=4: ref.append((float(r[0]), float(r[3]), float(r[4]) if len(r)>4 else 0.0))

# Match each REF point to the MC bin containing its xval; compute agreement.
rel=[]; chi2=0.0; n=0
print(f"  bin x      MC density   ALEPH ref    rel.dev")
for (xv, yref, yerr) in ref:
    md=None
    for (xc,xlo,xhi,dens) in mc:
        if xlo<=xv<xhi: md=dens; break
    if md is None or yref<=0: continue
    rd=(md-yref)/yref; rel.append(abs(rd)); n+=1
    if yerr>0: chi2+=((md-yref)/yerr)**2
    if n<=8 or abs(rd)>0.5:
        print(f"  {xv:.3f}    {md:.4e}   {yref:.4e}   {rd:+.1%}")
meanrel=sum(rel)/len(rel) if rel else float('nan')
print(f"\n  observable {suf}: {n} bins compared")
print(f"  mean |MC-ALEPH|/ALEPH = {meanrel:.1%}")
if n>0: print(f"  chi2/ndf (ALEPH errors only) = {chi2/n:.2f}")
