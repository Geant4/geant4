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
// $Id: EnergyHists.cc,v 1.1 2003-05-27 13:44:48 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "EnergyHists.hh"
#include "RunAction.hh"
#include "EnergyAngleCrossSection.hh"


EnergyHists::EnergyHists()
{}


EnergyHists::EnergyHists(G4double lo, G4double hi, G4double size)
  : loBinEdge(lo), hiBinEdge(hi), binSize(size)
{}


EnergyHists::~EnergyHists()
{}


void EnergyHists::CreateHists(IHistogramFactory* hFactory)
{
  G4int nBins = (G4int)((hiBinEdge - loBinEdge)/binSize);

  hists.push_back(hFactory->
    createHistogram1D("5 deg MC",  nBins, loBinEdge, hiBinEdge) );
  hists.push_back(hFactory->
    createHistogram1D("11 deg MC", nBins, loBinEdge, hiBinEdge) );
  hists.push_back(hFactory->
    createHistogram1D("15 deg MC", nBins, loBinEdge, hiBinEdge) );
  hists.push_back(hFactory->
    createHistogram1D("20 deg MC", nBins, loBinEdge, hiBinEdge) );
  hists.push_back(hFactory->
    createHistogram1D("30 deg MC", nBins, loBinEdge, hiBinEdge) );

  for (G4int i = 0; i < (G4int)hists.size(); i++) 
                                  hists[i]->fill(1000.0,1.0e-30);  
}


void EnergyHists::FillHists(G4double charge, G4double ke, 
                              G4double weight, G4int spectrum)
{
  hists[spectrum]->fill(ke,weight);
} 


void EnergyHists::ScaleHists(G4double factor, G4int spectrum)
{
  hists[spectrum]->scale(factor);
} 









