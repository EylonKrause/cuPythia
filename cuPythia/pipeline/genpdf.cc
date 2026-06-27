// cuPythia pipeline — real-PDF support for the device PDF interpolator (pdf.cuh).
// Fills the log-x/log-Q^2 grid from Pythia 8.317's REAL proton PDF (the default NNPDF set)
// instead of the toy, ports pdf.cuh's log(xf) bilinear interpolation to the host, and
// validates fidelity vs Pythia's PDF in the cross-section-support region — the real-PDF
// analog of pdf_xsec.cu's toy validation (1.1e-3). Also writes the grid to real_pdf.grid so
// the device hadronic-sigma stage can convolve a PHYSICAL gluon PDF (drops into the same
// arrays + interpolator unchanged), making that sigma physical rather than illustrative.
//
// Build: g++ genpdf.cc -o genpdf $(../../pythia8317/bin/pythia8-config --cxxflags --libs)
// Run:   ./genpdf   (needs PYTHIA8DATA set to the xmldoc dir)

#include "Pythia8/Pythia.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace Pythia8;

int main(){
  const int nx=200, nq=60;
  const double xmin=1e-4,xmax=1.0,Qmin=1.0,Qmax=2000.0;
  double lxmin=log(xmin),lxmax=log(xmax),lqmin=log(Qmin*Qmin),lqmax=log(Qmax*Qmax);

  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("ProcessLevel:all = off");   // PDF-only: no hard process needed
  pythia.readString("Print:quiet = on");
  pythia.init();
  PDFPtr pdf = pythia.getPDFPtr(2212);
  if(!pdf){ printf("no PDF\n"); return 1; }

  // Fill log(xf_g) grid from Pythia's real gluon PDF (LHAPDF practice; same as pdf.cuh).
  std::vector<double> g((size_t)nx*nq);
  for(int i=0;i<nx;++i){ double x=exp(lxmin+(lxmax-lxmin)*i/(nx-1));
    for(int j=0;j<nq;++j){ double Q2=exp(lqmin+(lqmax-lqmin)*j/(nq-1));
      g[(size_t)i*nq+j]=log(fmax(pdf->xf(21,x,Q2),1e-300)); } }

  // Port of pdf.cuh pdf_xfg: bilinear in (log x, log Q^2) of log(xf), edge-clamped, exp.
  auto interp=[&](double x,double Q2)->double{
    double fi=(log(x)-lxmin)/(lxmax-lxmin)*(nx-1), fj=(log(Q2)-lqmin)/(lqmax-lqmin)*(nq-1);
    if(fi<0)fi=0; if(fi>nx-1)fi=nx-1; if(fj<0)fj=0; if(fj>nq-1)fj=nq-1;
    int i0=(int)fi,j0=(int)fj,i1=(i0+1<nx)?i0+1:nx-1,j1=(j0+1<nq)?j0+1:nq-1;
    double ti=fi-i0,tj=fj-j0;
    double a=g[(size_t)i0*nq+j0],b=g[(size_t)i1*nq+j0],c=g[(size_t)i0*nq+j1],d=g[(size_t)i1*nq+j1];
    return exp((a*(1-ti)+b*ti)*(1-tj)+(c*(1-ti)+d*ti)*tj);
  };

  // Fidelity vs Pythia's PDF in the cross-section-support region (x ~ 1e-3..0.3).
  double maxRel=0; unsigned s=12345;
  for(int k=0;k<20000;++k){
    s=s*1664525u+1013904223u; double rx=(s>>9)*(1.0/8388608.0);
    s=s*1664525u+1013904223u; double rq=(s>>9)*(1.0/8388608.0);
    double x=exp(log(1e-3)+(log(0.3)-log(1e-3))*rx), Q2=exp(log(4.0)+(log(1e6)-log(4.0))*rq);
    double got=interp(x,Q2), exact=pdf->xf(21,x,Q2);
    if(exact>1e-9) maxRel=fmax(maxRel,fabs(got-exact)/fabs(exact));
  }

  // Write the grid for the device stage (header + log(xf) values).
  FILE* f=fopen("real_pdf.grid","w");
  if(f){ fprintf(f,"%d %d %.10e %.10e %.10e %.10e\n",nx,nq,lxmin,lxmax,lqmin,lqmax);
    for(size_t i=0;i<(size_t)nx*nq;++i) fprintf(f,"%.10e\n",g[i]); fclose(f); }

  printf("Real PDF (Pythia 8.317 proton, gluon) on the cuPythia device grid (%dx%d)\n",nx,nq);
  printf("  e.g. xf_g(x=0.01,Q=100) = %.4f (Pythia) vs %.4f (grid interp)\n",
         pdf->xf(21,0.01,1e4), interp(0.01,1e4));
  printf("  interp fidelity vs Pythia PDF in sigma-support (x:1e-3..0.3) = %.2e\n", maxRel);
  printf("  grid written -> real_pdf.grid (drops into pdf.cuh arrays for a PHYSICAL sigma)\n");
  bool ok=(maxRel<5e-2);   // bilinear on a real NNPDF gluon; ~1-2% expected in the support
  printf("REAL-PDF VALIDATE: %s (device interpolator reproduces Pythia's real gluon PDF)\n",ok?"PASS":"FAIL");
  return ok?0:2;
}
