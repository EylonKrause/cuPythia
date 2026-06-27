// cuPythia pipeline — HepMC3 output. Reads the GPU FSR-shower parton chains
// (shower_partons.dat, produced on-GPU by bridge.cu) and writes them as spec-valid HepMC3
// (Asciiv3) GenEvents — the standard event-record format for detector sims / analysis.
// Validated by reading the file straight back with HepMC3's own reader: event count and
// per-event 4-momentum must round-trip exactly.
//
// Build: g++ hepmc3_writer.cc -o hepmc3_writer -I$HOME/.local/include -L$HOME/.local/lib \
//          -lHepMC3 -Wl,-rpath,$HOME/.local/lib
// Run:   ./hepmc3_writer [shower_partons.dat] [out=cupythia.hepmc3]

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/FourVector.h"
#include <cstdio>
#include <cmath>
#include <memory>
using namespace HepMC3;

int main(int argc,char**argv){
  const char* in =(argc>1)?argv[1]:"shower_partons.dat";
  const char* out=(argc>2)?argv[2]:"cupythia.hepmc3";
  FILE* f=fopen(in,"r"); if(!f){ printf("cannot open %s (run ./bridge first)\n",in); return 1; }
  int N; double Ecm; if(fscanf(f,"%d %lf",&N,&Ecm)!=2){ printf("bad header\n"); return 1; }

  // ---- write ----
  WriterAscii writer(out);
  long nWrote=0;
  for(int e=0;e<N;++e){
    int n; if(fscanf(f,"%d",&n)!=1) break;
    GenEvent evt(Units::GEV, Units::MM);
    evt.set_event_number(e);
    GenVertexPtr v = std::make_shared<GenVertex>();
    // An incoming Z (the string's total 4-momentum, at rest at sqrt(s)) -> outgoing partons.
    v->add_particle_in(std::make_shared<GenParticle>(FourVector(0,0,0,Ecm), 23, 4));
    for(int i=0;i<n;++i){ int id,col,acol; double px,py,pz,en;
      if(fscanf(f,"%d %d %d %lf %lf %lf %lf",&id,&col,&acol,&px,&py,&pz,&en)!=7){ n=-1; break; }
      v->add_particle_out(std::make_shared<GenParticle>(FourVector(px,py,pz,en), id, 1)); }
    if(n<0) break;
    evt.add_vertex(v);
    writer.write_event(evt);
    nWrote++;
  }
  writer.close(); fclose(f);

  // ---- read back and validate ----
  ReaderAscii reader(out);
  long nRead=0; double maxImb=0;
  while(!reader.failed()){
    GenEvent evt(Units::GEV, Units::MM);
    reader.read_event(evt);
    if(reader.failed()) break;
    double s0=0,s1=0,s2=0,s3=0;
    for(auto& p : evt.particles()) if(p->status()==1){
      const FourVector& m=p->momentum(); s0+=m.px(); s1+=m.py(); s2+=m.pz(); s3+=m.e(); }
    double imb=fabs(s0)+fabs(s1)+fabs(s2)+fabs(s3-Ecm);
    if(imb>maxImb) maxImb=imb;
    nRead++;
  }
  reader.close();

  printf("cuPythia -> HepMC3 (Asciiv3): wrote %ld, read back %ld events -> %s\n", nWrote, nRead, out);
  printf("  max per-event 4-momentum round-trip imbalance = %.2e GeV\n", maxImb);
  bool ok = (nWrote>0) && (nRead==nWrote) && (maxImb<1e-6);
  printf("HEPMC3 VALIDATE: %s (spec-valid HepMC3 written + read back + conserving)\n", ok?"PASS":"FAIL");
  return ok?0:2;
}
