// cuPythia pipeline — device-resident PDF evaluator (build step 3).
//
// A parton distribution f(x,Q^2) is read from a (log x, log Q^2) grid in device
// memory by bilinear interpolation, with edge CLAMPING (low-x / low-Q^2 freezing)
// — the GAPS-v2 lesson: omit freezing and you get an unphysical boundary dip.
// Here the grid is filled from a physically-shaped TOY gluon PDF so the machinery
// is validatable analytically; a real LHAPDF .dat grid plugs into the SAME arrays
// and the SAME interpolator with no kernel change.
#pragma once
#include <cmath>
#include <vector>
#include <cuda_runtime.h>

// Toy gluon xf_g(x,Q^2) = x*f_g: low-x rise ~ x^-0.3, valence-like (1-x)^5, mild
// logarithmic Q^2 scaling. Illustrative shape only (NOT a real PDF fit).
__host__ __device__ inline double toy_xfg(double x,double Q2){
  double Q02=1.0; double q=1.0+0.3*log(Q2/Q02); if(q<0.2) q=0.2;
  return 4.0*pow(x,-0.3)*pow(1.0-x,5.0)*q;
}

struct PdfGrid { int nx,nq; double lxmin,lxmax,lqmin,lqmax; double* xf; }; // xf[ix*nq+iq]

__host__ inline PdfGrid buildToyGrid(int nx,int nq,double xmin,double xmax,double Qmin,double Qmax){
  PdfGrid g; g.nx=nx; g.nq=nq;
  g.lxmin=log(xmin); g.lxmax=log(xmax); g.lqmin=log(Qmin*Qmin); g.lqmax=log(Qmax*Qmax);
  std::vector<double> h((size_t)nx*nq);
  for(int i=0;i<nx;++i){ double x=exp(g.lxmin+(g.lxmax-g.lxmin)*i/(nx-1));
    for(int j=0;j<nq;++j){ double Q2=exp(g.lqmin+(g.lqmax-g.lqmin)*j/(nq-1));
      h[(size_t)i*nq+j]=log(fmax(toy_xfg(x,Q2),1e-300)); } }  // store log(xf) (LHAPDF practice)
  cudaMalloc(&g.xf,(size_t)nx*nq*sizeof(double));
  cudaMemcpy(g.xf,h.data(),(size_t)nx*nq*8,cudaMemcpyHostToDevice);
  return g;
}
__host__ inline void freeGrid(PdfGrid& g){ cudaFree(g.xf); }

// xf_g(x,Q^2) by bilinear interpolation in (log x, log Q^2), CLAMPED at edges (freezing).
__host__ __device__ inline double pdf_xfg(const PdfGrid& g,double x,double Q2){
  double fi=(log(x)-g.lxmin)/(g.lxmax-g.lxmin)*(g.nx-1);
  double fj=(log(Q2)-g.lqmin)/(g.lqmax-g.lqmin)*(g.nq-1);
  if(fi<0)fi=0; if(fi>g.nx-1)fi=g.nx-1;            // low/high-x freezing
  if(fj<0)fj=0; if(fj>g.nq-1)fj=g.nq-1;            // low/high-Q^2 freezing
  int i0=(int)fi, j0=(int)fj;
  int i1=(i0+1<g.nx)?i0+1:g.nx-1, j1=(j0+1<g.nq)?j0+1:g.nq-1;
  double ti=fi-i0, tj=fj-j0;
  double a=g.xf[(size_t)i0*g.nq+j0], b=g.xf[(size_t)i1*g.nq+j0];   // stored as log(xf)
  double c=g.xf[(size_t)i0*g.nq+j1], d=g.xf[(size_t)i1*g.nq+j1];
  double lxf=(a*(1-ti)+b*ti)*(1-tj)+(c*(1-ti)+d*ti)*tj;
  return exp(lxf);                                                 // interpolate in log(xf)
}
// f_g(x,Q^2) = xf/x.
__host__ __device__ inline double pdf_g(const PdfGrid& g,double x,double Q2){ return pdf_xfg(g,x,Q2)/x; }
