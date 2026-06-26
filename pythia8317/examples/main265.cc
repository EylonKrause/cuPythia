// main265.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>
//          Stephen Mrenna <mrenna@fnal.gov>
//          Philip Ilten <philten@cern.ch>

// Keywords: hadronization; reweighting; tuning

// Demonstrates how to reweight an event for different hadronization
// kinematics parameters, after the event has been produced, i.e. post-hoc
// rather than in-situ reweighting. The result here should be equivalent to
// the in-situ kinematics reweighitng of main263. The power of this method
// is that by saving a minimal set of information per event (10
// doubles per string), the entire event can be reweighted to whatever
// flavor parameters are requested by the user.

// Pythia includes.
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/KinematicVariations.h"
#include "Pythia8/PythiaParallel.h"

using namespace Pythia8;

int main() {

  // Define the new set of kinematic parameters that we wish to
  // reweight to. This can be set at the event level, but we define it
  // here so that we can compare with the in-situ reweighting.
  double aLund   = 0.6;  // StringZ:aLund, default 0.68
  double bLund   = 0.9;  // StringZ:bLund, default 0.98
  double rFactC  = 1.3;  // StringZ:rFactC, default 1.32
  double rFactB  = 0.9;  // StringZ:rFactB, default 0.855
  double ptSigma = 0.3;  // StringPT:sigma, default 0.335

  // Create, configure, and initialize Pythia.
  PythiaParallel pythia;
  pythia.readString("Beams:idA = 11");
  pythia.readString("Beams:idB = -11");
  pythia.readString("Beams:eCM = 91.189");
  pythia.readString("PDF:lepton = off");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");
  pythia.readString("HadronLevel:Decay = off");
  pythia.readString("StringFlav:popcornRate = 0");
  //pythia.readString("Parallelism:numThreads = 1");
  pythia.readString("Parallelism:processAsync = on");
  pythia.readString("VariationFrag:List = {kineVar frag:aLund="
    + to_string(aLund) + " frag:bLund=" + to_string(bLund) + " frag:rFactC="
    + to_string(rFactC) + " frag:rFactB=" + to_string(rFactB)
    + " frag:ptSigma=" + to_string(ptSigma) + "}");

  // Unlike flavor reweighting, the post-hoc reweighting must be
  // turned on explicitly, even if in-situ reweighting has been
  // enabled. This is because the information that must be stored for
  // kinematic reweighting is substantial.
  pythia.readString("VariationFrag:z = on");
  pythia.readString("VariationFrag:pT = on");
  pythia.init();

  // Define the plot title.
  string title = "default: (" +
      toString(pythia.settings.parm("StringZ:aLund")) + ", " +
      toString(pythia.settings.parm("StringZ:bLund")) + ", " +
      toString(pythia.settings.parm("StringZ:rFactC")) + ", " +
      toString(pythia.settings.parm("StringZ:rFactB")) + ", " +
      toString(pythia.settings.parm("StringPT:sigma")) + "), " +
      "variation: (" + toString(aLund) + ", " + toString(bLund) + ", " +
      toString(rFactC) + ", " + toString(rFactB) + ", " +
      toString(ptSigma) + ") ";

  // Create the reweighting tool given a settings object. Here, since
  // the parallel framework is being used, the settings from the
  // helper Pythia object which is used to intialize the other Pythia
  // versions is passed.
  map<int, KinematicVariations*> vars;
  pythia.foreach([&](Pythia* pythiaPtr) {
    vars[pythiaPtr->settings.mode("Parallelism:index")] =
      new KinematicVariations(pythiaPtr->settings);
  });
  // Alternatively, the tool can be created by passing the default
  // parameters.
  // KinematicVariations vars(aLundDefault, bLundDefault,
  //   rFactCDefault, rFactBDefault, ptSigmaDefault);

  // Define multiplicity histograms.
  // default:     the default parameters in Pythia
  // posthoc-wgt: post-hoc reweighted using the Pythia event
  // posthoc-str: post-hoc reweighted using the saved string break
  // insitu:      in-situ reweighted
  // rerun:       a run with the varied parameters
  vector<string> names = {
    "default", "posthoc-wgt", "posthoc-str", "insitu", "rerun"};
  map<string, Hist> hists;
  for (string &name : names)
    hists[name] = Hist(name, 25, 2, 51);

  // Track the weights.
  map<string, double> wgts, sumWgts, sumWgt2s;
  for (string &name : names)
    wgts[name] = sumWgts[name] = sumWgt2s[name] = 0;
  names.pop_back();

  // This mutual exclusion (mutex) object controls access to weights
  // and histograms.
  mutex wgtsMutex;

  // Run events.
  int nEvent = 1e6;
  // This defines a lambda function that acts as a callback.
  // This function is called for each event generated.
  // The argument is a pointer to the instance that generated the event.
  pythia.run(nEvent, [&](Pythia* pythiaPtr) {

    Event &event = pythiaPtr->event;
    int idx = pythiaPtr->settings.mode("Parallelism:index");

    // The necessary information for reweighting later can be saved to
    // a string. Note, other serialization implementations can be
    // used, and could then be implemented with symmetric
    // KinematicVariations::read and KinematicVariations::write methods.
    string zIntSaved = vars[idx]->write(pythiaPtr->info.weightContainerPtr
      ->weightsFragmentation.zIntBreaks);
    string zDblSaved = vars[idx]->write(pythiaPtr->info.weightContainerPtr
      ->weightsFragmentation.zDblBreaks);
    string pTSaved = vars[idx]->write(pythiaPtr->info.weightContainerPtr
      ->weightsFragmentation.pTBreaks);

    // Calculate the weight for the event, assuming we already have
    // the weight container and associated string breaks.
    double wgtPH = vars[idx]->weight(aLund, bLund, rFactC, rFactB, ptSigma,
      pythiaPtr->info.weightContainerPtr->weightsFragmentation.zIntBreaks,
      pythiaPtr->info.weightContainerPtr->weightsFragmentation.zDblBreaks,
      pythiaPtr->info.weightContainerPtr->weightsFragmentation.pTBreaks);

    // If instead we have saved the breaks to a string, as we did
    // above, we can calculate the weight from the saved string.
    double wgtPHStr = vars[idx]->weight(aLund, bLund, rFactC, rFactB, ptSigma,
      intVectorAttributeValue(zIntSaved, ""),
      doubleVectorAttributeValue(zDblSaved, ""),
      doubleVectorAttributeValue(pTSaved, ""));

    // We can also use the in-situ reweighting.
    double wgtIS = pythiaPtr->info.weightValueByIndex(
      pythiaPtr->info.numberOfWeights() - 1);

    // Lock for weights and histograms.
    wgtsMutex.lock();

    // For the default parameters, the weight is just 1.
    wgts["default"] = 1;
    wgts["posthoc-wgt"] = wgtPH;
    wgts["posthoc-str"] = wgtPHStr;
    wgts["insitu"] = wgtIS;

    // Keep track of the weights.
    for (string &name : names) {
      sumWgts[name]  += wgts[name];
      sumWgt2s[name] += pow2(wgts[name]);
    }

    // Fill the histograms.
    int mult = 0;
    for (const Particle &prt : event) {
      if (!prt.isFinal()) continue;
      if (prt.isCharged()) ++mult;
    }
    for (string &name : names) hists[name].fill(mult, wgts[name]);
    wgtsMutex.unlock();
  });
  pythia.stat();

  // Rerun Pythia with the varied parameters.
  pythia.settings.parm("StringZ:aLund",  aLund);
  pythia.settings.parm("StringZ:bLund",  bLund);
  pythia.settings.parm("StringZ:rFactC", rFactC);
  pythia.settings.parm("StringZ:rFactB", rFactB);
  pythia.settings.parm("StringPT:sigma", ptSigma);
  pythia.settings.wvec("VariationFrag:List", {});
  pythia.settings.flag("VariationFrag:z", false);
  pythia.settings.flag("VariationFrag:pT", false);

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  pythia.run( nEvent, [&](Pythia* pythiaPtr) {

    Event &event = pythiaPtr->event;
    int mult = 0;
    for (const Particle &prt : event) {
      if (!prt.isFinal()) continue;
      if (prt.isCharged()) ++mult;
    }
    wgtsMutex.lock();
    sumWgts["rerun"]  += 1;
    sumWgt2s["rerun"] += 1;
    hists["rerun"].fill(mult, 1.);
    wgtsMutex.unlock();
  });
  pythia.stat();

  // Normalize the histograms.
  for (auto &hist : hists) hist.second /= sumWgts[hist.first];

  // Print the histogram ratios.
  string xlabel = "multiplicity";
  for (auto &hist : hists)
    cout << "\n" << hist.first << hist.second/hists["default"];

  // Print the reweighting stats.
  // The 1 - mu should be statistically consistent with zero if the
  // reweighting has proper coverage.
  // The n_eff gives the statistical power of the reweighted sample.
  for (string &name : names) {
    double w(sumWgts[name]), w2(sumWgt2s[name]), n(sumWgts["default"]);
    cout << name << "\n"
         << "\t1 - mu = " << scientific << setprecision(3) << abs(1. - w/n)
         << " +- "<< abs(1. - sqrt((w2/n - pow2(w/n))*n/(n - 1)))/sqrt(n)
         << "\n\tn_eff  = " << scientific << setprecision(3) << w/sqrt(n*w2)
         << "\n";
  }

  // Create Python plot.
  HistPlot hpl("plot265");
  hpl.frame("fig265", title, xlabel, "n(variation)/n(default)");
  for (string &name : names)
    hpl.add(hists[name]/hists["default"], "e", name);
  hpl.add(hists["rerun"]/hists["default"], "e", "rerun");
  hpl.plot();

  // Clean up the memory.
  for (auto var : vars) {
    delete var.second;
  }

  return 0;
}
