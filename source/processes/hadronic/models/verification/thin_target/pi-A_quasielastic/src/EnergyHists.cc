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
//
// $Id: EnergyHists.cc,v 1.1 2003-07-31 01:21:09 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "EnergyHists.hh"
#include "RunAction.hh"
#include "EnergyAngleCrossSection.hh"

#include "G4ios.hh"


EnergyHists::EnergyHists()
{}


EnergyHists::EnergyHists(G4double lo, G4double hi, G4double size)
  : loBinEdge(lo), hiBinEdge(hi), binSize(size)
{}


EnergyHists::~EnergyHists()
{}


void EnergyHists::CreateHists()
{
  G4int nBins = (G4int)((hiBinEdge - loBinEdge)/binSize);

  G4int i;
  for (i = 0; i < nBins; i++) hist19_5deg.push_back(0);
  hists.push_back(&hist19_5deg);
  for (i = 0; i < nBins; i++) hist30deg.push_back(0);
  hists.push_back(&hist30deg);
  for (i = 0; i < nBins; i++) hist42deg.push_back(0);
  hists.push_back(&hist42deg);
}


void EnergyHists::FillHists(G4double charge, G4double ke, G4int spectrum)
{
  if (charge < 0) {
    G4int nBins = (G4int)((hiBinEdge - loBinEdge)/binSize);
    G4int chan = (G4int)(ke/binSize);
    if (chan < nBins) {
      G4double binContents = (*hists[spectrum])[chan];
      binContents += 1.0;
      (*hists[spectrum])[chan] = binContents;
    }
  }
} 


void EnergyHists::ScaleHists(G4double factor, G4int spectrum)
{
  G4int nBins = (G4int)((hiBinEdge - loBinEdge)/binSize);

  G4String fileName; 
  if (spectrum == 0) {
    fileName = "Pi-quasi_821MeV_19_5deg.hist" ;
  } else if (spectrum == 1) {
    fileName = "Pi-quasi_821MeV_30deg.hist" ;
  } else if (spectrum == 2) {
    fileName = "Pi-quasi_821MeV_42deg.hist" ;
  }
  ofstream fcout(fileName, ios::app);
  fcout << " Pion KE bin center (MeV)   cross section(mb/sr/MeV)  stat. error (mb/sr/MeV)" << endl;

  for (G4int chan = 0; chan < nBins; chan++) {
    G4double binContents = (*hists[spectrum])[chan];
    G4double csValue = binContents*factor;
    G4double csError = sqrt(binContents)*factor;
    G4int binCenter = (G4int)(binSize*(0.5 + (G4int)chan));

    fcout << "               " << binCenter 
          << "                      " << csValue 
          << "                      " << csError  << endl;
  }
} 




