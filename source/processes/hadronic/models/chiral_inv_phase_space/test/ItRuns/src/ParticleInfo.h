#ifndef ANAParticleInfo_h
#define ANAParticleInfo_h
#include <fstream>
#include "g4std/vector"
#include "Analysis/src/Particle.h"
#include "Analysis/src/Plot.h"
#include "Analysis/src/DataPoint.h"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4KaonPlus.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4KaonMinus.hh"
#include "G4Triton.hh"
#include "globals.hh"
#include "Analysis/src/VFilter.h"

class ANAParticleInfo
{
  public:
    ANAParticleInfo(G4double xSec, G4String aFileName);
    ~ANAParticleInfo() {} // Needs to clean up memory.
    void Analyse();
    
  private:
    G4std::vector<ANAParticle> theParicles;
    // pdg, info
    // info = lowBin, xsec/totalXsec
    vector<VANAPlot *> thePlots;
};

ANAParticleInfo::ANAParticleInfo(G4double xSec, G4String aFileName)
{
  ifstream theData(aFileName);
  for(;;)
  {
    ANAParticle aPart;
    if(!aPart.Init(theData)) break;
    theParicles.push_back(aPart);
  }
    
  ANAPlot<ANADataPoint> * aNewPlot = 0;
// ANAPlot<ANADataPoint> ::ANAPlot<DataPoint> (G4int aPDG, G4double aMass, G4double aTotalXsec, G4String fn)
// pi+
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(211, G4PionPlus::PionPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piplus"));
  thePlots.push_back(aNewPlot);
// pi-
  aNewPlot = new ANAPlot<ANADataPoint>(-211, G4PionMinus::PionMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piminus"));
  thePlots.push_back(aNewPlot);
// p 
  aNewPlot = new ANAPlot<ANADataPoint>(2212, G4Proton::ProtonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".proton"));
  thePlots.push_back(aNewPlot);
// deuteron
  aNewPlot = new ANAPlot<ANADataPoint>(1002, G4Deuteron::DeuteronDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".deuteron")); // 1000*Z+A
  thePlots.push_back(aNewPlot);
// K+
  aNewPlot = new ANAPlot<ANADataPoint>(321, G4KaonPlus::KaonPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonplus"));
  thePlots.push_back(aNewPlot);
// triton
  aNewPlot = new ANAPlot<ANADataPoint>(1003, G4Triton::TritonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".triton")); // 1000*Z+A
  thePlots.push_back(aNewPlot);
// he3
  aNewPlot = new ANAPlot<ANADataPoint>(2003, G4He3::He3Definition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he3")); // 1000*Z+A
  thePlots.push_back(aNewPlot);
// he4
  aNewPlot = new ANAPlot<ANADataPoint>(2004, G4Alpha::AlphaDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he4")); // 1000*Z+A
  thePlots.push_back(aNewPlot);
// K-
  aNewPlot = new ANAPlot<ANADataPoint>(-321, G4KaonMinus::KaonMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonminus"));
  thePlots.push_back(aNewPlot);
}

inline
void ANAParticleInfo::Analyse()
{
  G4int aParticle;
  G4int aPlot;
  for(aParticle=0; aParticle<theParicles.size(); aParticle++)
  {
    G4int pdg = theParicles[aParticle].GetPDGCode();
    G4double energy = theParicles[aParticle].GetEnergy();
    G4double weight = theParicles[aParticle].GetWeight();
    for(aPlot = 0; aPlot<thePlots.size(); aPlot++)
    {
      if(thePlots[aPlot]->Insert(pdg, energy, weight)) break;
    }
  }
  for(aPlot = 0; aPlot<thePlots.size(); aPlot++)
  {
    G4cout << "New plot:"<<G4endl;
    thePlots[aPlot]->DumpInfo(G4cout);
  }
}

#endif
