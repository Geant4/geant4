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

//#define debug

#include <fstream>
#include "g4std/vector"
#include "GammaNuclear/src/Particle.h"
#include "GammaNuclear/src/Plot.h"
#include "GammaNuclear/src/DataPoint.h"
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
#include "GammaNuclear/src/VFilter.h"
#include "GammaNuclear/src/DoubleBandFilter.h"

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
  G4double cBin;
  G4double raGr = 3.14159265/180. ; 
  G4double hBin = raGr*2. ; 

  cBin  = raGr*57.;
  DoubleBandFilter * bin1 = new DoubleBandFilter(cos(cBin+hBin), cos(cBin-hBin), ".57degrees");

  cBin  = raGr*77.;
  DoubleBandFilter * bin2 = new DoubleBandFilter(cos(cBin+hBin), cos(cBin-hBin), ".77degrees");

  cBin  = raGr*97.;
  DoubleBandFilter * bin3 = new DoubleBandFilter(cos(cBin+hBin), cos(cBin-hBin), ".97degrees");

  cBin  = raGr*117.;
  DoubleBandFilter * bin4 = new DoubleBandFilter(cos(cBin+hBin), cos(cBin-hBin), ".117degrees");

  cBin  = raGr*127.;
  DoubleBandFilter * bin5 = new DoubleBandFilter(cos(cBin+hBin), cos(cBin-hBin), ".127degrees");

  ANAPlot<ANADataPoint, TVANAFilter<G4double> > * aNewPlot = 0;

// ANAPlot<ANADataPoint> ::ANAPlot<DataPoint> (G4int aPDG, G4double aMass, G4double aTotalXsec, G4String fn)
// protons -------------------- 
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
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2212, G4Proton::ProtonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".proton"), bin5);
  thePlots.push_back(aNewPlot);
// deuteron ------------------
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
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1002, G4Deuteron::DeuteronDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".deuteron"), bin5); // 1000*Z+A
  thePlots.push_back(aNewPlot);
// triton --------------------
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
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(1003, G4Triton::TritonDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".triton"), bin5); // 1000*Z+A
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
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2003, G4He3::He3Definition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he3"), bin5); // 1000*Z+A
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
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2004, G4Alpha::AlphaDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".he4"), bin5); // 1000*Z+A
  thePlots.push_back(aNewPlot);
}

inline
void ANAParticleInfo::Analyse()
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
    ProcessOne(aPart);
  }
}

void ANAParticleInfo::Plot(G4String aPreFix, G4int aStatistics)
{
  for(int aPlot = 0; aPlot<thePlots.size(); aPlot++)
  {
#ifdef debug
    G4cout << "New plot:"<<G4endl;
#endif
    thePlots[aPlot]->SetNevents(aStatistics);
    thePlots[aPlot]->DumpInfo(G4cout, aPreFix);
  }
}

inline void ANAParticleInfo::ProcessOne(ANAParticle aPart)
{
    G4int pdg = aPart.GetPDGCode();
    G4double energy = aPart.GetEnergy();
    G4double weight = aPart.GetWeight();
    for(int aPlot = 0; aPlot<thePlots.size(); aPlot++)
    {
      if(thePlots[aPlot]->Filter( &(aPart) ))
      {
        if(thePlots[aPlot]->Insert(pdg, energy, weight)) break;
      }
    }  
}
#endif
