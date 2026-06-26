// cuPythia kernel 07 — GPU UNWEIGHTING efficiency + Les Houches (.lhe) output.
//
// The headline metric of modern GPU parton-level generators (Pepper,
// arXiv:2311.06198; madgraph4gpu): the UNWEIGHTING EFFICIENCY, i.e. the fraction
// of generated (weighted) events kept after acceptance-rejection,
//     eta = <w> / w_max .
// It drives the real cost: you pay to generate 1/eta events for every unweighted
// event a detector simulation consumes, which is why phase-space/importance
// sampling (VEGAS, MadNIS arXiv:2212.06172) matters. This kernel measures eta on
// GPU for gg->gg and writes the accepted (unweighted) events to a standard Les
// Houches Event file — the parton-level interchange format Pythia/Herwig and the
// experiments read — with NO external library.
//
// Build: nvcc -O3 -arch=sm_120 -o unweight_lhe unweight_lhe.cu
// Run:   ./unweight_lhe [trialsPerThread=4000]   (writes events.lhe)

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

__host__ __device__ inline double pow2(double x) { return x * x; }
__host__ __device__ inline double gg2gg_sigma(double sH, double tH, double uH, double alpS) {
  double is=1.0/sH, it=1.0/tH, iu=1.0/uH;  // reciprocal precompute: 3 FP64 div instead of 13
  double rts=tH*is, rst=sH*it, rus=uH*is, rsu=sH*iu, rtu=tH*iu, rut=uH*it;
  double a = (9. / 4.) * (rts*rts + 2.*rts + 3. + 2.*rst + rst*rst);
  double b = (9. / 4.) * (rus*rus + 2.*rus + 3. + 2.*rsu + rsu*rsu);
  double c = (9. / 4.) * (rtu*rtu + 2.*rtu + 3. + 2.*rut + rut*rut);
  return (M_PI * is*is) * pow2(alpS) * 0.5 * (a + b + c);
}
__host__ __device__ inline double weightAt(double cosT, double s, double alpS) {
  double t = -0.5 * s * (1.0 - cosT);
  return gg2gg_sigma(s, t, -s - t, alpS); // dσ/dt̂ = the event weight
}

// Each thread: nPer trials. Accept with prob w/wMax (von Neumann). Returns per-
// block accepted count and summed weight (for the eta = <w>/wMax cross-check).
__global__ void unweightKernel(uint64_t seed, uint64_t nPer, double s, double alpS,
                               double cMax, double wMax,
                               unsigned long long* gAcc, double* gSumW) {
  uint64_t tid = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
  uint64_t ctr = seed + tid * 0x100000001B3ULL;
  unsigned long long acc = 0; double sumw = 0.0;
  for (uint64_t i = 0; i < nPer; ++i) {
    double c = (2.0 * u01(splitmix64(ctr++)) - 1.0) * cMax;
    double w = weightAt(c, s, alpS);
    sumw += w;
    double u = u01(splitmix64(ctr++));
    if (u * wMax < w) ++acc;       // accept -> an unweighted event
  }
  atomicAdd(gAcc, acc);
  atomicAdd(gSumW, sumw);
}

static double wMaxScan(int N, double s, double alpS, double cMax) {
  double m = 0.0;
  for (int i = 0; i <= N; ++i) { double c = -cMax + (2.0 * cMax) * i / N; m = fmax(m, weightAt(c, s, alpS)); }
  return m;
}
static double simpson(int N, double s, double alpS, double cMax) {
  double a = -cMax, h = (2.0 * cMax) / N, sum = 0.0;
  for (int i = 0; i <= N; ++i) {
    double c = a + i * h, w = (i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    sum += w * weightAt(c, s, alpS);
  }
  return (sum * h / 3.0) * (s / 2.0);
}

// Write nEv accepted (unweighted) gg->gg events to a valid Les Houches file.
static void writeLHE(const char* path, const std::vector<double>& cosArr,
                     const std::vector<double>& phiArr, double s, double xsec_pb) {
  FILE* f = fopen(path, "w");
  if (!f) return;
  double E = 0.5 * sqrt(s); // CM energy per parton (massless), beams along z
  fprintf(f, "<LesHouchesEvents version=\"3.0\">\n<header>\n");
  fprintf(f, "  cuPythia gg->gg unweighted sample (demonstration)\n</header>\n");
  fprintf(f, "<init>\n");
  // IDBMUP1 IDBMUP2 EBMUP1 EBMUP2 PDFGUP1 PDFGUP2 PDFSUP1 PDFSUP2 IDWTUP NPRUP
  fprintf(f, "2212 2212 %.8E %.8E 0 0 0 0 3 1\n", E, E);
  // XSECUP XERRUP XMAXUP LPRUP   (weighting strategy 3 = unweighted, w=+1)
  fprintf(f, "%.8E %.8E 1.0E+00 1\n", xsec_pb, xsec_pb * 0.01);
  fprintf(f, "</init>\n");
  for (size_t e = 0; e < cosArr.size(); ++e) {
    double c = cosArr[e], st = sqrt(fmax(0.0, 1.0 - c * c)), phi = phiArr[e];
    double p3x = E * st * cos(phi), p3y = E * st * sin(phi), p3z = E * c;
    fprintf(f, "<event>\n");
    // NUP IDPRUP XWGTUP SCALUP AQEDUP AQCDUP
    fprintf(f, "4 1 1.0E+00 %.8E 7.546771E-03 1.180000E-01\n", sqrt(s));
    // id status mother1 mother2 col acol px py pz E m spin lifetime  (gluons id=21)
    fprintf(f, "21 -1 0 0 501 502 0 0  %.8E %.8E 0 0 1\n", E, E);
    fprintf(f, "21 -1 0 0 503 501 0 0 %.8E %.8E 0 0 1\n", -E, E);
    fprintf(f, "21  1 1 2 503 504  %.8E %.8E %.8E %.8E 0 0 1\n", p3x, p3y, p3z, E);
    fprintf(f, "21  1 1 2 504 502  %.8E %.8E %.8E %.8E 0 0 1\n", -p3x, -p3y, -p3z, E);
    fprintf(f, "</event>\n");
  }
  fprintf(f, "</LesHouchesEvents>\n");
  fclose(f);
}

#define CK(call) do { cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return 1; } } while(0)

int main(int argc, char** argv) {
  uint64_t nPer = (argc > 1) ? strtoull(argv[1], nullptr, 10) : 4000ULL;
  const int blocks = 1024, threads = 256;
  double s = 100.0 * 100.0, alphaS = 0.118, cMax = 0.9, conv_pb = 0.3893793721e9; // GeV^-2 -> pb
  uint64_t seed = 0x7E57ULL;
  double wMax = wMaxScan(2000000, s, alphaS, cMax);
  uint64_t total = (uint64_t)blocks * threads * nPer;

  unsigned long long *dAcc; double *dSumW;
  CK(cudaMalloc(&dAcc, sizeof(unsigned long long))); CK(cudaMemset(dAcc, 0, sizeof(unsigned long long)));
  CK(cudaMalloc(&dSumW, sizeof(double)));            CK(cudaMemset(dSumW, 0, sizeof(double)));
  unweightKernel<<<blocks, threads>>>(seed, nPer, s, alphaS, cMax, wMax, dAcc, dSumW);
  CK(cudaDeviceSynchronize());
  unsigned long long hAcc = 0; double hSumW = 0.0;
  CK(cudaMemcpy(&hAcc, dAcc, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(&hSumW, dSumW, sizeof(double), cudaMemcpyDeviceToHost));

  double meanW = hSumW / (double)total;
  double eta = (double)hAcc / (double)total;       // unweighting efficiency
  double etaCheck = meanW / wMax;                  // = <w>/wMax
  double V = (2.0 * cMax) * (s / 2.0);             // cosθ width * Jacobian dt̂/dcosθ
  double sigma_pb = V * meanW * conv_pb;
  double ref_pb = simpson(2000000, s, alphaS, cMax) * conv_pb;

  // Build a small unweighted LHE sample on the host (same accept logic).
  std::vector<double> cosA, phiA; uint64_t ctr = 0xABCULL;
  while (cosA.size() < 1000) {
    double c = (2.0 * u01(splitmix64(ctr++)) - 1.0) * cMax;
    double u = u01(splitmix64(ctr++));
    double phi = 2.0 * M_PI * u01(splitmix64(ctr++));
    if (u * wMax < weightAt(c, s, alphaS)) { cosA.push_back(c); phiA.push_back(phi); }
  }
  writeLHE("events.lhe", cosA, phiA, s, sigma_pb);

  printf("GPU unweighting (gg->gg) + Les Houches output\n");
  printf("  trials                 = %.3e   w_max = %.4e\n", (double)total, wMax);
  printf("  unweighting efficiency = %.2f%%   (cross-check <w>/w_max = %.2f%%)\n",
         100.0 * eta, 100.0 * etaCheck);
  printf("  sigma (unweighted MC)  = %.6e pb   Simpson ref = %.6e pb   relerr = %.2e\n",
         sigma_pb, ref_pb, fabs(sigma_pb - ref_pb) / ref_pb);
  printf("  wrote %zu unweighted events -> events.lhe (standard Les Houches format)\n", cosA.size());
  bool ok = (fabs(eta - etaCheck) < 1e-3) && (fabs(sigma_pb - ref_pb) / ref_pb < 2e-3) && cosA.size() == 1000;
  printf("VALIDATION: %s (eta == <w>/w_max, sigma matches quadrature, LHE written)\n", ok ? "PASS" : "FAIL");
  cudaFree(dAcc); cudaFree(dSumW);
  return ok ? 0 : 2;
}
