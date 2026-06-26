// cuPythia — device-resident event record (Structure-of-Arrays). THE DATA PLANE.
//
// A real GPU generator (vs the isolated demo kernels) keeps an event on the device
// across every stage — hard process -> PDF -> reweight -> shower -> hadronization ->
// decays — reading/writing particles HERE. Per-event randomness is counter-based
// (rng.cuh) so any single event is O(1)-reproducible on any node, which makes the
// bit-identical CPU-equivalence test (the validation cornerstone, per GAPS) trivial.
//
// Layout: particle p of event e is at flat index e*maxPart + p (SoA -> coalesced
// when one thread owns one event). Fields chosen per the architecture research
// (GAPS arXiv:2403.08692, MadtRex arXiv:2510.05100, MCnet weight naming 2203.08230).
//
// Status (Pythia-like, simplified): -21 incoming beam parton; 23 outgoing hard
// parton; 51 shower-emitted parton; 1 final-state; 2 decayed/intermediate.
#pragma once
#include <cstdint>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

struct DeviceEvents {
  int nEvents, maxPart, nVar;            // nVar = weight variations (>=1; slot 0 = nominal)
  // --- per particle [nEvents*maxPart] ---
  double  *px,*py,*pz,*e,*m;             // four-momentum + mass
  int     *pdg,*status,*col,*acol,*mo1,*mo2,*d1,*d2;   // id, status, colour, mother+daughter links
  // --- per event [nEvents] ---
  int     *nPart;                        // live particle count
  uint64_t*seed;                         // RNG key
  double  *scale,*x1,*x2;                // hard scale; incoming parton momentum fractions
  int     *flavA,*flavB;                 // incoming parton flavours
  unsigned char *active;                 // GAPS active/endShower flag (1 = still evolving)
  double  *weight;                       // [nEvents*nVar] weight vector (k=0 nominal)
};

__host__ inline DeviceEvents allocEvents(int nEvents, int maxPart, int nVar=1) {
  DeviceEvents ev; ev.nEvents=nEvents; ev.maxPart=maxPart; ev.nVar=nVar;
  size_t np=(size_t)nEvents*maxPart, ne=(size_t)nEvents;
  auto D=[](size_t bytes){ void* p=nullptr; cudaMalloc(&p,bytes); return p; };
  ev.px=(double*)D(np*8); ev.py=(double*)D(np*8); ev.pz=(double*)D(np*8);
  ev.e =(double*)D(np*8); ev.m =(double*)D(np*8);
  ev.pdg=(int*)D(np*4); ev.status=(int*)D(np*4); ev.col=(int*)D(np*4); ev.acol=(int*)D(np*4);
  ev.mo1=(int*)D(np*4); ev.mo2=(int*)D(np*4); ev.d1=(int*)D(np*4); ev.d2=(int*)D(np*4);
  ev.nPart=(int*)D(ne*4); ev.seed=(uint64_t*)D(ne*8);
  ev.scale=(double*)D(ne*8); ev.x1=(double*)D(ne*8); ev.x2=(double*)D(ne*8);
  ev.flavA=(int*)D(ne*4); ev.flavB=(int*)D(ne*4);
  ev.active=(unsigned char*)D(ne);
  ev.weight=(double*)D(ne*(size_t)nVar*8);
  cudaMemset(ev.nPart,0,ne*4);
  cudaMemset(ev.active,1,ne);            // all events start active
  return ev;
}
__host__ inline void freeEvents(DeviceEvents& ev){
  cudaFree(ev.px);cudaFree(ev.py);cudaFree(ev.pz);cudaFree(ev.e);cudaFree(ev.m);
  cudaFree(ev.pdg);cudaFree(ev.status);cudaFree(ev.col);cudaFree(ev.acol);
  cudaFree(ev.mo1);cudaFree(ev.mo2);cudaFree(ev.d1);cudaFree(ev.d2);
  cudaFree(ev.nPart);cudaFree(ev.seed);cudaFree(ev.scale);cudaFree(ev.x1);cudaFree(ev.x2);
  cudaFree(ev.flavA);cudaFree(ev.flavB);cudaFree(ev.active);cudaFree(ev.weight);
}

// Append a particle to event e. Caller owns event e (one thread/event) -> no atomics.
// Returns the new slot index, or -1 if the event is full.
__device__ inline int addParticle(DeviceEvents& ev,int e,
                                   double px,double py,double pz,double en,double mass,
                                   int pdg,int status,int col,int acol,
                                   int mo1,int mo2,int d1=-1,int d2=-1){
  int p=ev.nPart[e]; if(p>=ev.maxPart) return -1;
  size_t i=(size_t)e*ev.maxPart+p;
  ev.px[i]=px; ev.py[i]=py; ev.pz[i]=pz; ev.e[i]=en; ev.m[i]=mass;
  ev.pdg[i]=pdg; ev.status[i]=status; ev.col[i]=col; ev.acol[i]=acol;
  ev.mo1[i]=mo1; ev.mo2[i]=mo2; ev.d1[i]=d1; ev.d2[i]=d2;
  ev.nPart[e]=p+1; return p;
}
