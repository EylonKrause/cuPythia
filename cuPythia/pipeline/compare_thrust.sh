#!/bin/bash
# Overlay cuPythia GPU FSR thrust vs Pythia SimpleTimeShower thrust, with per-bin ratio
# and a shape-agreement metric (mean |ratio-1| over populated bins, and chi2-like).
cd "$(dirname "$0")"
G=thrust_gpu.dat; P="${1:-thrust_pythia.dat}"
[ -f "$G" ] || { echo "missing $G (run ./shower_fsr)"; exit 1; }
[ -f "$P" ] || { echo "missing $P (run ./thrust_pythia)"; exit 1; }
echo "  GPU FSR shower  vs  $P"
echo "  (1-T) bin        GPU density    Pythia density   ratio G/P"
echo "  ----------------------------------------------------------"
paste <(grep -v '^#' "$G") <(grep -v '^#' "$P") | awk '
{
  lo=$1; hi=$2; g=$3; p=$6;
  if (p>0 && g>0) { r=g/p; n++; sumdev+=(r>1?r-1:1-r); chi+=(g-p)*(g-p)/p; }
  else r=0;
  printf "  %.3f-%.3f      %10.4e     %10.4e     %s\n", lo, hi, g, p, (r>0?sprintf("%.3f",r):"-");
}
END {
  printf "\n  populated bins: %d\n", n;
  printf "  mean |ratio-1| over populated bins = %.1f%%\n", 100*sumdev/n;
  printf "  (lower = better shape agreement)\n";
}'
