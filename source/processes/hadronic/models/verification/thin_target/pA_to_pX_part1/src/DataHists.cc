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
// $Id: DataHists.cc,v 1.1 2003-05-27 13:44:48 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "DataHists.hh"
#include "EnergyAngleCrossSection.hh"


DataHists::DataHists()
{}


DataHists::DataHists(G4double lo, G4double hi, G4double size)
  : loBinEdge(lo), hiBinEdge(hi), binSize(size)
{}


DataHists::~DataHists()
{}


void DataHists::CreateHists(IHistogramFactory* hFactory)
{
  G4int nBins = (G4int)((hiBinEdge - loBinEdge)/binSize);

  datahists.push_back(hFactory->
           createHistogram1D("5 deg data",  nBins, loBinEdge, hiBinEdge) );
  datahists.push_back(hFactory->
           createHistogram1D("11 deg data", nBins, loBinEdge, hiBinEdge) );
  datahists.push_back(hFactory->
           createHistogram1D("15 deg data", nBins, loBinEdge, hiBinEdge) );
  datahists.push_back(hFactory->
           createHistogram1D("20 deg data", nBins, loBinEdge, hiBinEdge) );
  datahists.push_back(hFactory->
           createHistogram1D("30 deg data", nBins, loBinEdge, hiBinEdge) );
}


void DataHists::PlotData(EnergyAngleCrossSection* theData)
{
  G4int Nhists = datahists.size();
  for (G4int k = 0; k < Nhists; k++) datahists[k]->fill(1000.0, 1.0e-20);

  ///////////////////////////////////////
  // (p,p') cross sections (mb/sr/MeV)
  ///////////////////////////////////////

  G4std::vector<G4double>& angles = theData->GetAngles();
  G4int Nangles = angles.size();

  if (Nangles != Nhists) 
     G4Exception(" Number of angles does not match number of histograms ");

  for (G4int i = 0; i < Nangles; i++) {
    G4std::vector<G4std::pair<G4double, G4double> >& keSpectrum =
                            theData->GetEnergySpectrum(angles[i]);

    // Loop over energies
    //
    for (G4int j = 0; j < (G4int)keSpectrum.size(); j++) {
      datahists[i]->fill(keSpectrum[j].first,0.001*keSpectrum[j].second);
    }
  }
} 









