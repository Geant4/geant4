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
#include "Analysis/src/DoubleBandFilter.h"

class ANAParticleInfo
{
  public:
    ANAParticleInfo(G4double xSec, G4String aFileName);
    ~ANAParticleInfo() {} // Needs to clean up memory.
    void Analyse(G4String aPreFix);
    
  private:
    // pdg, info
    // info = lowBin, xsec/totalXsec
    vector<VANAPlot *> thePlots;
    G4String theFileName;
};

ANAParticleInfo::ANAParticleInfo(G4double xSec, G4String aFileName) : theFileName(aFileName)
{
   
  G4double aMean;
  G4double halfBin = 0.07; 

//  DoubleBandFilter * bin1 = new DoubleBandFilter(0.25,  0.4, ".70degree"); // cos(th)
  aMean = cos(3.14159265*70./180.);
  DoubleBandFilter * bin1 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".70degree"); // cos(th)

//  DoubleBandFilter * bin2 = new DoubleBandFilter(-0.25, 0.25, ".90degree"); // cos(th)
  aMean = cos(3.14159265*90./180.);
  DoubleBandFilter * bin2 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".90degree"); // cos(th)

//  DoubleBandFilter * bin3 = new DoubleBandFilter(-0.7, -0.25, ".118degree"); // cos(th)
  aMean = cos(3.14159265*118./180.);
  DoubleBandFilter * bin3 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".118degree"); // cos(th)

//  DoubleBandFilter * bin4 = new DoubleBandFilter(-1.0,  -0.7, ".160degree"); // cos(th)
  aMean = cos(3.14159265*160./180.);
  DoubleBandFilter * bin4 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".160degree"); // cos(th)
 
  ANAPlot<ANADataPoint, TVANAFilter<G4double> > * aNewPlot = 0;
// ANAPlot<ANADataPoint> ::ANAPlot<DataPoint> (G4int aPDG, G4double aMass, G4double aTotalXsec, G4String fn)
// pi+
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(211, G4PionPlus::PionPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piplus"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(211, G4PionPlus::PionPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piplus"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(211, G4PionPlus::PionPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piplus"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(211, G4PionPlus::PionPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piplus"), bin4);
  thePlots.push_back(aNewPlot);
// pi-
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(-211, G4PionMinus::PionMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piminus"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(-211, G4PionMinus::PionMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piminus"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(-211, G4PionMinus::PionMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piminus"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(-211, G4PionMinus::PionMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".piminus"), bin4);
  thePlots.push_back(aNewPlot);
// p 
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2212, G4Proton::ProtonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".proton"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2212, G4Proton::ProtonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".proton"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2212, G4Proton::ProtonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".proton"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2212, G4Proton::ProtonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".proton"), bin4);
  thePlots.push_back(aNewPlot);
// deuteron
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1002, G4Deuteron::DeuteronDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".deuteron"), bin1); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1002, G4Deuteron::DeuteronDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".deuteron"), bin2); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1002, G4Deuteron::DeuteronDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".deuteron"), bin3); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1002, G4Deuteron::DeuteronDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".deuteron"), bin4); // 1000*Z+A
  thePlots.push_back(aNewPlot);
// K+
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(321, G4KaonPlus::KaonPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonplus"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(321, G4KaonPlus::KaonPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonplus"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(321, G4KaonPlus::KaonPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonplus"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(321, G4KaonPlus::KaonPlusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonplus"), bin4);
  thePlots.push_back(aNewPlot);
// triton
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1003, G4Triton::TritonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".triton"), bin1); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1003, G4Triton::TritonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".triton"), bin2); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1003, G4Triton::TritonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".triton"), bin3); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1003, G4Triton::TritonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".triton"), bin4); // 1000*Z+A
  thePlots.push_back(aNewPlot);
// he3
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2003, G4He3::He3Definition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he3"), bin1); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2003, G4He3::He3Definition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he3"), bin2); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2003, G4He3::He3Definition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he3"), bin3); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2003, G4He3::He3Definition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he3"), bin4); // 1000*Z+A
  thePlots.push_back(aNewPlot);
// he4
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2004, G4Alpha::AlphaDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he4"), bin1); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2004, G4Alpha::AlphaDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he4"), bin2); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2004, G4Alpha::AlphaDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he4"), bin3); // 1000*Z+A
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2004, G4Alpha::AlphaDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he4"), bin4); // 1000*Z+A
  thePlots.push_back(aNewPlot);
// K-
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(-321, G4KaonMinus::KaonMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonminus"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(-321, G4KaonMinus::KaonMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonminus"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(-321, G4KaonMinus::KaonMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonminus"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(-321, G4KaonMinus::KaonMinusDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".kaonminus"), bin4);
  thePlots.push_back(aNewPlot);
}

inline
void ANAParticleInfo::Analyse(G4String aPreFix)
{
  G4int aParticle;
  G4int aPlot;
  ifstream theData(theFileName);
  G4int counter = 0;
  for(;;)
  {
    counter++;
    if(counter == 10000*(counter/10000)) 
       G4cout << "taken care of "<<counter<<" particles." << G4endl;
    ANAParticle aPart;
    if(!aPart.Init(theData)) break;
    G4int pdg = aPart.GetPDGCode();
    G4double energy = aPart.GetEnergy();
    G4double weight = aPart.GetWeight();
    for(aPlot = 0; aPlot<thePlots.size(); aPlot++)
    {
      if(thePlots[aPlot]->Filter( &(aPart) ))
      {
        if(thePlots[aPlot]->Insert(pdg, energy, weight)) break;
      }
    }
  }
  for(aPlot = 0; aPlot<thePlots.size(); aPlot++)
  {
    G4cout << "New plot:"<<G4endl;
    thePlots[aPlot]->DumpInfo(G4cout, aPreFix);
  }
}

#endif
