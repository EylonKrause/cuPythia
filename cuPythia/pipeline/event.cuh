// cuPythia — device-resident event record (Structure-of-Arrays). THE DATA PLANE.
//
// The whole point of a real GPU generator (vs the isolated demo kernels) is that
// an event never leaves the device between stages. Every pipeline stage — hard
// process -> parton shower -> hadronization -> decays — reads and writes its
// particles HERE, in device memory. Per-event randomness is counter-based
// (rng.cuh) so any single event is O(1)-reproducible on any node.
//
// Layout: particle p of event e is at flat index  e*maxPart + p  (SoA: when one
// thread owns one event, same-field accesses across events are coalesced).
//
// Status (Pythia-like, simplified): -21 incoming beam parton; 23 outgoing
// hard-process parton; 1 final-state particle; 2 intermediate/decayed.
#pragma once
#include <cstdint>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

struct DeviceEvents {
  int nEvents, maxPart;
  double  *px,*py,*pz,*e,*m;              // four-momentum + mass  [nEvents*maxPart]
  int     *pdg,*status,*col,*acol,*mo1,*mo2;
  int     *nPart;                         // live particle count per event [nEvents]
  uint64_t*seed;                          // per-event RNG key      [nEvents]
  double  *weight,*scale;                 // per-event weight, hard scale [nEvents]
};

__host__ inline DeviceEvents allocEvents(int nEvents, int maxPart) {
  DeviceEvents ev; ev.nEvents = nEvents; ev.maxPart = maxPart;
  size_t np = (size_t)nEvents * maxPart;
  auto D = [](size_t bytes){ void* p=nullptr; cudaMalloc(&p, bytes); return p; };
  ev.px=(double*)D(np*8); ev.py=(double*)D(np*8); ev.pz=(double*)D(np*8);
  ev.e =(double*)D(np*8); ev.m =(double*)D(np*8);
  ev.pdg=(int*)D(np*4); ev.status=(int*)D(np*4); ev.col=(int*)D(np*4);
  ev.acol=(int*)D(np*4); ev.mo1=(int*)D(np*4); ev.mo2=(int*)D(np*4);
  ev.nPart=(int*)D((size_t)nEvents*4);
  ev.seed=(uint64_t*)D((size_t)nEvents*8);
  ev.weight=(double*)D((size_t)nEvents*8); ev.scale=(double*)D((size_t)nEvents*8);
  cudaMemset(ev.nPart, 0, (size_t)nEvents*4);
  return ev;
}
__host__ inline void freeEvents(DeviceEvents& ev) {
  cudaFree(ev.px); cudaFree(ev.py); cudaFree(ev.pz); cudaFree(ev.e); cudaFree(ev.m);
  cudaFree(ev.pdg); cudaFree(ev.status); cudaFree(ev.col); cudaFree(ev.acol);
  cudaFree(ev.mo1); cudaFree(ev.mo2); cudaFree(ev.nPart); cudaFree(ev.seed);
  cudaFree(ev.weight); cudaFree(ev.scale);
}

// Append a particle to event e. Caller owns event e (one thread per event) so no
// atomics are needed; returns the new particle's slot, or -1 if the event is full.
__device__ inline int addParticle(DeviceEvents& ev, int e,
                                   double px,double py,double pz,double en,double mass,
                                   int pdg,int status,int col,int acol,int mo1,int mo2) {
  int p = ev.nPart[e];
  if (p >= ev.maxPart) return -1;
  size_t i = (size_t)e * ev.maxPart + p;
  ev.px[i]=px; ev.py[i]=py; ev.pz[i]=pz; ev.e[i]=en; ev.m[i]=mass;
  ev.pdg[i]=pdg; ev.status[i]=status; ev.col[i]=col; ev.acol[i]=acol;
  ev.mo1[i]=mo1; ev.mo2[i]=mo2;
  ev.nPart[e] = p + 1;
  return p;
}
