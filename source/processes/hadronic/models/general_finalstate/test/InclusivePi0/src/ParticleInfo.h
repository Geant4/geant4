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
#include "InclusivePi0/src/Particle.h"
#include "InclusivePi0/src/Plot.h"
#include "InclusivePi0/src/DataPoint.h"
#include "G4PionZero.hh"
#include "globals.hh"
#include "InclusivePi0/src/VFilter.h"
#include "InclusivePi0/src/DoubleBandFilter.h"

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
  G4double halfBin = 0.05; 

  aMean = .15;
  DoubleBandFilter * bin1 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".xl.15"); 
  aMean = .25;
  DoubleBandFilter * bin2 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".xl.25"); 
  aMean = .35;
  DoubleBandFilter * bin3 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".xl.35"); 
  aMean = .45;
  DoubleBandFilter * bin4 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".xl.45"); 
  aMean = .55;
  DoubleBandFilter * bin5 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".xl.55"); 
  aMean = .65;
  DoubleBandFilter * bin6 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".xl.65"); 
  aMean = .75;
  DoubleBandFilter * bin7 = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".xl.75"); 
  ANAPlot<ANADataPoint, TVANAFilter<G4double> > * aNewPlot = 0;
// ANAPlot<ANADataPoint> ::ANAPlot<DataPoint> (G4int aPDG, G4double aMass, G4double aTotalXsec, G4String fn)
// pi+
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(111, G4PionZero::PionZeroDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".pizero"), bin1);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(111, G4PionZero::PionZeroDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".pizero"), bin2);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(111, G4PionZero::PionZeroDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".pizero"), bin3);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(111, G4PionZero::PionZeroDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".pizero"), bin4);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(111, G4PionZero::PionZeroDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".pizero"), bin5);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(111, G4PionZero::PionZeroDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".pizero"), bin6);
  thePlots.push_back(aNewPlot);
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(111, G4PionZero::PionZeroDefinition()->GetPDGMass(), 
                                       xSec, aFileName+G4String(".pizero"), bin7);
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
    for(int aPlot = 0; aPlot<thePlots.size(); aPlot++)
    {
      if(thePlots[aPlot]->Filter( &(aPart) ))
      {
        if(thePlots[aPlot]->Insert(pdg, energy, weight)) break;
      }
    }  
}
#endif
