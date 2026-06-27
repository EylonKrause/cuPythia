#!/usr/bin/env python3
# Validate the baryon mechanism from a cuPythia hadron dump (hadronize_mr ... out.dat): per-event
# electric-charge and baryon-number conservation (must be exactly 0 for an e+e- -> q qbar event),
# and the inclusive baryon rates (p+pbar, Lambda+Lbar per event) vs the LEP/PDG world averages.
# 4-momentum conservation is already checked in-kernel; this is the FLAVOUR-side correctness gate.
import sys
# 3*charge of the particle named by the POSITIVE pdg code (antiparticle negates); baryon if |code|>1000.
Q3 = {211:3,111:0,321:3,311:0,130:0,310:0,221:0,331:0,113:0,213:3,223:0,333:0,313:0,323:3,22:0,
      # charm/bottom mesons (from g->qqbar endpoints): D,D*,Ds,Ds*,B,B*,Bs,Bs*,Bc,quarkonia
      411:3,421:0,431:3,413:3,423:0,433:3,441:0,443:0,
      511:0,521:3,531:0,513:0,523:3,533:0,541:3,543:3,551:0,553:0,
      2212:3,2112:0,3122:0,3222:3,3212:0,3112:-3,3312:-3,3322:0,3334:-3,
      1114:-3,2114:0,2214:3,2224:6,3114:-3,3214:0,3224:3,3314:-3,3324:0}
UNK=set()
def q3(pid):
  a=abs(pid); s=1 if pid>0 else -1
  if a not in Q3: UNK.add(a)
  return s*Q3.get(a,0)
def bnum(pid): a=abs(pid); return (1 if pid>0 else -1) if a>1000 else 0

f=open(sys.argv[1]); N,Ecm=f.readline().split(); N=int(N)
badQ=badB=0; nEvt=0; nBar=0; nHad=0
cnt={}
for _ in range(N):
    line=f.readline()
    if not line: break
    n=int(line); sq=sb=0
    for _ in range(n):
        p=f.readline().split(); pid=int(p[0])
        sq+=q3(pid); sb+=bnum(pid)
        if abs(pid)>1000: nBar+=1
        cnt[abs(pid)]=cnt.get(abs(pid),0)+1
    nHad+=n; nEvt+=1
    if sq!=0: badQ+=1
    if sb!=0: badB+=1
pbar=cnt.get(2212,0)/nEvt; lam=cnt.get(3122,0)/nEvt
print(f"events {nEvt}, mean hadrons {nHad/nEvt:.2f}, baryons/evt {nBar/nEvt:.3f} ({100.0*nBar/nHad:.1f}% of hadrons)")
print(f"  charge non-conserving events : {badQ}/{nEvt}")
print(f"  baryon-number non-conserving : {badB}/{nEvt}")
print(f"  p+pbar/evt = {pbar:.3f}  (PDG LEP ~1.05)   Lambda+Lbar/evt = {lam:.3f}  (PDG LEP ~0.39)")
if UNK: print(f"  WARNING unknown ids (charge assumed 0): {sorted(UNK)}")
ok = (badQ==0 and badB==0)
print(f"BARYON-FLAVOUR VALIDATION: {'PASS' if ok else 'FAIL'} (charge + baryon number conserved every event)")
sys.exit(0 if ok else 2)
