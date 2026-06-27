// Validate cuPythia's 2-loop running alpha_s (the -DAS_2LOOP path in shower_inc.cuh) against
// Pythia 8.317's own AlphaStrong order=2, and EXTRACT the Lambda_{3,4,5} constants to hardcode.
//
// Method: replicate Pythia's AlphaStrong::init Lambda-finding (StandardModel.cc:66-140, NITER=10)
// and alphaS evaluation (StandardModel.cc:189-247) in pure math, then compare alpha_s(mu^2) at a
// spread of scales to Pythia's AlphaStrong (same inputs: alpha_s(MZ)=0.1365, mc=1.5, mb=4.8).
//
// Build: g++ -O2 -std=c++17 as2loop_validate.cc -o as2loop_validate \
//   -I<pythia>/include -L<pythia>/lib -Wl,-rpath,<pythia>/lib -lpythia8 -ldl
// Run:   ./as2loop_validate
#include "Pythia8/StandardModel.h"
#include <cstdio>
#include <cmath>
using namespace Pythia8;

static const double MZ=91.188, MB=4.8, MC=1.5, ASMZ=0.1365;
static const double b15=348.0/529.0, b14=462.0/625.0, b13=64.0/81.0;

// Replicated Pythia 2-loop Lambda-finding (asMZ-fixed -> compute once).
static double L5g, L4g, L3g;
void findLambdas(){
  L5g=MZ*exp(-6.0*M_PI/(23.0*ASMZ));
  for(int i=0;i<10;++i){ double ls=2.0*log(MZ/L5g),lls=log(ls); double corr=1.0-b15*lls/ls;
    double vi=ASMZ/corr; L5g=MZ*exp(-6.0*M_PI/(23.0*vi)); }
  double lsB=2.0*log(MB/L5g),llsB=log(lsB); double vB=12.0*M_PI/(23.0*lsB)*(1.0-b15*llsB/lsB);
  L4g=L5g;
  for(int i=0;i<10;++i){ double ls=2.0*log(MB/L4g),lls=log(ls); double corr=1.0-b14*lls/ls;
    double vi=vB/corr; L4g=MB*exp(-6.0*M_PI/(25.0*vi)); }
  double lsC=2.0*log(MC/L4g),llsC=log(lsC); double vC=12.0*M_PI/(25.0*lsC)*(1.0-b14*llsC/lsC);
  L3g=L4g;
  for(int i=0;i<10;++i){ double ls=2.0*log(MC/L3g),lls=log(ls); double corr=1.0-b13*lls/ls;
    double vi=vC/corr; L3g=MC*exp(-6.0*M_PI/(27.0*vi)); }
}
// Replicated 2-loop evaluation.
double as2loop(double mu2){
  double s2min=1.33*1.33*L3g*L3g; if(mu2<s2min) mu2=s2min;  // Pythia SAFETYMARGIN2 freeze
  double Lam2,b0,b1;
  if(mu2>MB*MB){ Lam2=L5g*L5g; b0=23.0; b1=b15; }
  else if(mu2>MC*MC){ Lam2=L4g*L4g; b0=25.0; b1=b14; }
  else { Lam2=L3g*L3g; b0=27.0; b1=b13; }
  double ls=log(mu2/Lam2);
  return 12.0*M_PI/(b0*ls)*(1.0-b1*log(ls)/ls);
}

int main(){
  findLambdas();
  AlphaStrong as; as.setThresholds(MC,MB,171.0); as.init(ASMZ,2,5,false);
  printf("Lambda (mine vs Pythia):  L5 %.8f / %.8f   L4 %.8f / %.8f   L3 %.8f / %.8f\n",
         L5g,as.Lambda5(), L4g,as.Lambda4(), L3g,as.Lambda3());
  double scales[]={0.25,1.0,2.25,5.0,23.04,100.0,MZ*MZ/4.0,MZ*MZ,2000.0};
  double maxRel=0;
  printf("  scale2      mine        Pythia      relErr\n");
  for(double s2:scales){ double a=as2loop(s2), p=as.alphaS(s2); double r=fabs(a-p)/p;
    if(r>maxRel)maxRel=r; printf("  %9.3f  %.6f   %.6f   %.2e\n",s2,a,p,r); }
  double aMZ=as2loop(MZ*MZ);
  printf("alpha_s(M_Z^2) recovered = %.6f (target %.4f)\n", aMZ, ASMZ);
  bool ok=(maxRel<1e-4) && (fabs(aMZ-ASMZ)<1e-3);
  printf("VALIDATION: %s (2-loop alpha_s vs Pythia AlphaStrong order=2, max relErr %.2e)\n",
         ok?"PASS":"FAIL", maxRel);
  printf("HARDCODE: L5=%.8f L4=%.8f L3=%.8f\n", L5g,L4g,L3g);
  return ok?0:2;
}
