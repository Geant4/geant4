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
#include "g4std/vector"
#include "HadronNuclear/batch_monitor/Particle.h"
#include "HadronNuclear/batch_monitor/Plot.h"
#include "HadronNuclear/batch_monitor/DataPoint.h"
#include "G4Neutron.hh"
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
#include "HadronNuclear/batch_monitor/VFilter.h"
#include "HadronNuclear/batch_monitor/DoubleBandFilter.h"

class ANAParticleInfo
{
  public:
    ANAParticleInfo(G4double xSec, G4String aFileName);
    ~ANAParticleInfo() {} // Needs to clean up memory.
    void Analyse();
    void Plot(G4String aPreFix, G4int aStatistics);
    void ProcessOne(ANAParticle aPart);
    void makePlot(G4double bin, G4String degree, G4String file, G4double xSec);
    void makeProtonPlot(G4double bin, G4String degree, G4String file, G4double xSec);
    void makeInclusiveProtonPlot(G4String file, G4double xSec);
    void makeInclusiveNeutronPlot(G4String file, G4double xSec);
   
  private:
    // pdg, info
    // info = lowBin, xsec/totalXsec
    vector<VANAPlot *> thePlots;
    G4String theFileName;
};

ANAParticleInfo::ANAParticleInfo(G4double xSec, G4String aFileName) : theFileName(aFileName)
{
   
  makeInclusiveProtonPlot(aFileName+G4String(".proton"), xSec);
  makeInclusiveNeutronPlot(aFileName+G4String(".neutron"), xSec);

  makePlot(0.,   ".0degree",   aFileName+G4String(".neutron"), xSec);
  makePlot(7.5,  ".7.5degree", aFileName+G4String(".neutron"), xSec);
  makePlot(10.,  ".10degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(11.,  ".11degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(20.,  ".20degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(24.,  ".24degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(25.,  ".25degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(30.,  ".30degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(35.,  ".35degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(40.,  ".40degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(45.,  ".45degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(56.,  ".56degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(60.,  ".60degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(69.,  ".69degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(70.,  ".70degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(80.,  ".80degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(82.,  ".82degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(95.,  ".95degree",  aFileName+G4String(".neutron"), xSec);
  makePlot(100., ".100degree", aFileName+G4String(".neutron"), xSec);
  makePlot(106., ".106degree", aFileName+G4String(".neutron"), xSec);
  makePlot(120., ".120degree", aFileName+G4String(".neutron"), xSec);
  makePlot(121., ".121degree", aFileName+G4String(".neutron"), xSec);
  makePlot(133., ".133degree", aFileName+G4String(".neutron"), xSec);
  makePlot(134., ".134degree", aFileName+G4String(".neutron"), xSec);
  makePlot(140., ".140degree", aFileName+G4String(".neutron"), xSec);
  makePlot(144., ".144degree", aFileName+G4String(".neutron"), xSec);
  makePlot(145., ".145degree", aFileName+G4String(".neutron"), xSec);
  makePlot(150., ".150degree", aFileName+G4String(".neutron"), xSec);
  makePlot(160., ".160degree", aFileName+G4String(".neutron"), xSec);
  
  makeProtonPlot(24.,  ".24degree",  aFileName+G4String(".proton"), xSec);
  makeProtonPlot(25.,  ".25degree",  aFileName+G4String(".proton"), xSec);
  makeProtonPlot(35.,  ".35degree",  aFileName+G4String(".proton"), xSec);
  makeProtonPlot(45.,  ".45degree",  aFileName+G4String(".proton"), xSec);
  makeProtonPlot(56.,  ".56degree",  aFileName+G4String(".proton"), xSec);
  makeProtonPlot(69.,  ".69degree",  aFileName+G4String(".proton"), xSec);
  makeProtonPlot(70.,  ".70degree",  aFileName+G4String(".proton"), xSec);
  makeProtonPlot(82.,  ".82degree",  aFileName+G4String(".proton"), xSec);
  makeProtonPlot(95.,  ".95degree",  aFileName+G4String(".proton"), xSec);
  makeProtonPlot(106., ".106degree", aFileName+G4String(".proton"), xSec);
  makeProtonPlot(120., ".120degree", aFileName+G4String(".proton"), xSec);
  makeProtonPlot(121., ".121degree", aFileName+G4String(".proton"), xSec);
  makeProtonPlot(134., ".134degree", aFileName+G4String(".proton"), xSec);
  makeProtonPlot(140., ".140degree", aFileName+G4String(".proton"), xSec);
  makeProtonPlot(145., ".145degree", aFileName+G4String(".proton"), xSec);
  makeProtonPlot(150., ".150degree", aFileName+G4String(".proton"), xSec);

}

inline
void ANAParticleInfo::makePlot(double bin, G4String degree, G4String file, G4double xSec)
{
  G4double aMean = cos(bin);
  G4double halfBin = 0.07; 
  DoubleBandFilter * cut = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, degree); // cos(th)
  ANAPlot<ANADataPoint, TVANAFilter<G4double> > *aNewPlot;
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2112, G4Neutron::NeutronDefinition()->GetPDGMass(), 
                                       xSec, file, cut);
  thePlots.push_back(aNewPlot);
}

inline
void ANAParticleInfo::makeProtonPlot(double bin, G4String degree, G4String file, G4double xSec)
{
  G4double aMean = cos(bin);
  G4double halfBin = 0.07; 
  DoubleBandFilter * cut = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, degree); // cos(th)
  ANAPlot<ANADataPoint, TVANAFilter<G4double> > *aNewPlot;
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2212, G4Proton::ProtonDefinition()->GetPDGMass(), 
                                       xSec, file, cut);
  thePlots.push_back(aNewPlot);
}

inline
void ANAParticleInfo::makeInclusiveNeutronPlot(G4String file, G4double xSec)
{
  G4double aMean = 0;
  G4double halfBin = 10.; 
  DoubleBandFilter * cut = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".integral"); // cos(th)
  ANAPlot<ANADataPoint, TVANAFilter<G4double> > *aNewPlot;
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2112, G4Neutron::NeutronDefinition()->GetPDGMass(), 
                                       xSec, file, cut);
  thePlots.push_back(aNewPlot);
}

inline
void ANAParticleInfo::makeInclusiveProtonPlot(G4String file, G4double xSec)
{
  G4double aMean = 0;
  G4double halfBin = 10.; 
  DoubleBandFilter * cut = new DoubleBandFilter(aMean-halfBin, aMean+halfBin, ".integral"); // cos(th)
  ANAPlot<ANADataPoint, TVANAFilter<G4double> > *aNewPlot;
  aNewPlot = new ANAPlot<ANADataPoint, TVANAFilter<G4double> >(2212, G4Proton::ProtonDefinition()->GetPDGMass(), 
                                       xSec, file, cut);
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
        thePlots[aPlot]->Insert(pdg, energy, weight);
      }
    }  
}
#endif
