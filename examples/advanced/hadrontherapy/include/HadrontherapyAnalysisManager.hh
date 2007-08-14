//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// ----------------------------------------------------------------------------
// $Id: HadrontherapyAnalysisManager.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifdef G4ANALYSIS_USE
#ifndef HADRONTHERAPYANALYSISMANAGER_HH
#define HADRONTHERAPYANALYSISMANAGER_HH 1

#include "globals.hh"
# include <AIDA/AIDA.h>

namespace AIDA{
  class ITree; 
  class IAnalysisFactory;
  class ITreeFactory;
}

class HadrontherapyAnalysisManager
{
private:
  HadrontherapyAnalysisManager();

public:
  ~HadrontherapyAnalysisManager();
  
  static HadrontherapyAnalysisManager* getInstance();
  
  void book();
  // Book the histograms and ntuples in a .hbk file
  
  void FillEnergyDeposit(G4int voxelXId, G4int voxelYId, G4int voxelZId, 
                         G4double energyDeposit);
  // Fill the ntuple with the energy deposit in the phantom 

  void BraggPeak(G4int, G4double);
  // Fill 1D histogram with the Bragg peak in the phantom

  void SecondaryProtonEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary protons

   void SecondaryNeutronEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary neutrons

  void SecondaryAlphaEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary alpha particles

  void SecondaryGammaEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary gamma

  void SecondaryElectronEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary electrons

  void SecondaryTritonEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary tritons

  void SecondaryDeuteronEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary deuterons

  void SecondaryPionEnergyDeposit(G4int slice, G4double energy);
  // Fill 1D histogram with the energy deposit of secondary pions

  void electronEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary electrons originated in the phantom

  void gammaEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary gamma originated in the phantom

  void deuteronEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary deuterons originated in the phantom

  void tritonEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary tritons originated in the phantom

  void alphaEnergyDistribution(G4double secondaryParticleKineticEnergy);
  // Energy distribution of secondary alpha originated in the phantom

  void genericIonInformation(G4int, G4double, G4int, G4double);
 
  void finish();
  // Close the .hbk file with the histograms and the ntuples

private:
  static HadrontherapyAnalysisManager* instance;
  AIDA::IAnalysisFactory* aFact;
  AIDA::ITree* theTree; 
  AIDA::IHistogramFactory *histFact;
  AIDA::ITupleFactory *tupFact;
  AIDA::IHistogram1D *h1;
  AIDA::IHistogram1D *h2;
  AIDA::IHistogram1D *h3;
  AIDA::IHistogram1D *h4;
  AIDA::IHistogram1D *h5;
  AIDA::IHistogram1D *h6;
  AIDA::IHistogram1D *h7; 
  AIDA::IHistogram1D *h8; 
  AIDA::IHistogram1D *h9;
  AIDA::IHistogram1D *h10;
  AIDA::IHistogram1D *h11;
  AIDA::IHistogram1D *h12; 
  AIDA::IHistogram1D *h13; 
  AIDA::IHistogram1D *h14;
  AIDA::ITuple *ntuple;
  AIDA::ITuple *ionTuple;
};
#endif
#endif

