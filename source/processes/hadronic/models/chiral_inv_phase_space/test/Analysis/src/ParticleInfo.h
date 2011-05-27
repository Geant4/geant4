//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef ANAParticleInfo_h
#define ANAParticleInfo_h
#include <fstream>
#include <vector>
#include "Particle.h"
#include "Plot.h"
#include "DataPoint.h"
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
#include "VFilter.h"
#include "DoubleBandFilter.h"

using namespace std;

class ANAParticleInfo
{
  public:
    ANAParticleInfo(G4double xSec, G4String aFileName);
    ~ANAParticleInfo() {} // Needs to clean up memory.
    void Analyse();
    void Plot(G4String aPreFix, G4int aStatistics);
    void ProcessOne(ANAParticle aPart);
    
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

  aMean = cos(pi*70./180.);
  DoubleBandFilter* bin1 = new DoubleBandFilter(aMean-halfBin,aMean+halfBin,".70degree");

  aMean = cos(pi*90./180.);
  DoubleBandFilter* bin2 = new DoubleBandFilter(aMean-halfBin,aMean+halfBin,".90degree");

  aMean = cos(pi*118./180.);
  DoubleBandFilter* bin3 = new DoubleBandFilter(aMean-halfBin,aMean+halfBin,".118degree");

  aMean = cos(pi*137./180.);
  DoubleBandFilter* bin4 = new DoubleBandFilter(aMean-halfBin,aMean+halfBin,".137degree");
 
  aMean = cos(pi*160./180.);
  DoubleBandFilter* bin5 = new DoubleBandFilter(aMean-halfBin,aMean+halfBin,".160degree");

  ANAPlot<ANADataPoint, TVANAFilter<G4double> >* aNewPlot = 0;

// pi+
  G4double hMass=G4PionPlus::PionPlusDefinition()->GetPDGMass();
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 ( 211, hMass, xSec, aFileName+G4String(".piplus"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 ( 211, hMass, xSec, aFileName+G4String(".piplus"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 ( 211, hMass, xSec, aFileName+G4String(".piplus"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 ( 211, hMass, xSec, aFileName+G4String(".piplus"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 ( 211, hMass, xSec, aFileName+G4String(".piplus"), bin5);
  thePlots.push_back(aNewPlot);
// pi-
  hMass=G4PionMinus::PionMinusDefinition()->GetPDGMass();
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (-211, hMass, xSec, aFileName+G4String(".piminus"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (-211, hMass, xSec, aFileName+G4String(".piminus"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (-211, hMass, xSec, aFileName+G4String(".piminus"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (-211, hMass, xSec, aFileName+G4String(".piminus"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (-211, hMass, xSec, aFileName+G4String(".piminus"), bin5);
  thePlots.push_back(aNewPlot);
// p
  hMass=G4Proton::ProtonDefinition()->GetPDGMass();
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (2212, hMass, xSec, aFileName+G4String(".proton"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (2212, hMass, xSec, aFileName+G4String(".proton"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (2212, hMass, xSec, aFileName+G4String(".proton"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (2212, hMass, xSec, aFileName+G4String(".proton"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (2212, hMass, xSec, aFileName+G4String(".proton"), bin5);
  thePlots.push_back(aNewPlot);
// deuteron
  hMass=G4Deuteron::DeuteronDefinition()->GetPDGMass();
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                          (1000010020, hMass, xSec, aFileName+G4String(".deuteron"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                          (1000010020, hMass, xSec, aFileName+G4String(".deuteron"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                          (1000010020, hMass, xSec, aFileName+G4String(".deuteron"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                          (1000010020, hMass, xSec, aFileName+G4String(".deuteron"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                          (1000010020, hMass, xSec, aFileName+G4String(".deuteron"), bin5);
  thePlots.push_back(aNewPlot);
// K+
  hMass=G4KaonPlus::KaonPlusDefinition()->GetPDGMass();
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (321, hMass, xSec, aFileName+G4String(".kaonplus"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (321, hMass, xSec, aFileName+G4String(".kaonplus"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (321, hMass, xSec, aFileName+G4String(".kaonplus"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (321, hMass, xSec, aFileName+G4String(".kaonplus"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                                 (321, hMass, xSec, aFileName+G4String(".kaonplus"), bin5);
  thePlots.push_back(aNewPlot);
// triton
  hMass=G4Triton::TritonDefinition()->GetPDGMass();
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                            (1000010030, hMass, xSec, aFileName+G4String(".triton"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                            (1000010030, hMass, xSec, aFileName+G4String(".triton"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                            (1000010030, hMass, xSec, aFileName+G4String(".triton"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                            (1000010030, hMass, xSec, aFileName+G4String(".triton"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                            (1000010030, hMass, xSec, aFileName+G4String(".triton"), bin5);
  thePlots.push_back(aNewPlot);
// he3
  hMass=G4He3::He3Definition()->GetPDGMass();
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                              (10000020030, hMass, xSec, aFileName+G4String(".he3"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                              (10000020030, hMass, xSec, aFileName+G4String(".he3"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                              (10000020030, hMass, xSec, aFileName+G4String(".he3"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                              (10000020030, hMass, xSec, aFileName+G4String(".he3"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                              (10000020030, hMass, xSec, aFileName+G4String(".he3"), bin5);
  thePlots.push_back(aNewPlot);
// he4
  hMass=G4Alpha::AlphaDefinition()->GetPDGMass();
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (1000020040, hMass, xSec, aFileName+G4String(".he4"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (1000020040, hMass, xSec, aFileName+G4String(".he4"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (1000020040, hMass, xSec, aFileName+G4String(".he4"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (1000020040, hMass, xSec, aFileName+G4String(".he4"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (1000020040, hMass, xSec, aFileName+G4String(".he4"), bin5);
  thePlots.push_back(aNewPlot);
// K-
  hMass=G4KaonMinus::KaonMinusDefinition()->GetPDGMass();
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (-321, hMass, xSec, aFileName+G4String(".kaonminus"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (-321, hMass, xSec, aFileName+G4String(".kaonminus"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (-321, hMass, xSec, aFileName+G4String(".kaonminus"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (-321, hMass, xSec, aFileName+G4String(".kaonminus"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >
                               (-321, hMass, xSec, aFileName+G4String(".kaonminus"), bin5);
  thePlots.push_back(aNewPlot);
}

inline void ANAParticleInfo::Analyse()
{
  ifstream theData(theFileName);
  G4int counter = 0;
  for(;;)
  {
    counter++;
    if(!(counter%10000)) G4cout<<counter<<" particles are treated"<< G4endl;
    ANAParticle aPart;
    if(!aPart.Init(theData)) break;
    ProcessOne(aPart);
  }
}

inline void ANAParticleInfo::Plot(G4String aPreFix, G4int aStatistics)
{
  for(unsigned aPlot = 0; aPlot<thePlots.size(); ++aPlot)
  {
    G4cout << "New plot:"<<G4endl;
    thePlots[aPlot]->SetNevents(aStatistics);
    thePlots[aPlot]->DumpInfo(G4cout, aPreFix);
  }
}

inline void ANAParticleInfo::ProcessOne(ANAParticle aPart)
{
    G4int pdg = aPart.GetPDGCode();
    G4double energy = aPart.GetEnergy();
    G4double weight = aPart.GetWeight();
    for(unsigned aPlot = 0; aPlot<thePlots.size(); ++aPlot)
      if(thePlots[aPlot]->Filter( &(aPart) ) &&
         thePlots[aPlot]->Insert(pdg, energy, weight)) break;
}
#endif
