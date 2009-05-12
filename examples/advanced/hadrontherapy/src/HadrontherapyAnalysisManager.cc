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
// $Id: HadrontherapyAnalisysManager.cc;
// See more at: http://geant4infn.wikispaces.com/HadrontherapyExample

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

#ifdef  G4ANALYSIS_USE
#include "HadrontherapyAnalysisManager.hh"

HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::instance = 0;

/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisManager::HadrontherapyAnalysisManager() : 
  aFact(0), theTree(0), histFact(0), tupFact(0), h1(0), h2(0), h3(0),
  h4(0), h5(0), h6(0), h7(0), h8(0), h9(0), h10(0), h11(0), h12(0), h13(0), h14(0), ntuple(0),
  ionTuple(0)
{  
}
/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisManager::~HadrontherapyAnalysisManager() 
{ 
  delete ionTuple;
  ionTuple = 0;
  
  delete ntuple;
  ntuple = 0;

  delete h14;
  h14 = 0;

  delete h13;
  h13 = 0;

  delete h12;
  h12 = 0;

  delete h11;
  h11 = 0;

  delete h10;
  h10 = 0;

  delete h9;
  h9 = 0;

  delete h8;
  h8 = 0;

  delete h7;
  h7 = 0;

  delete h6;
  h6 = 0;

  delete h5;
  h5 = 0;

  delete h4;
  h4 = 0;

  delete h3;
  h3 = 0;

  delete h2;
  h2 = 0;

  delete h1;
  h1 = 0;
  
  delete tupFact;
  tupFact = 0;

  delete histFact;
  histFact = 0;

  delete theTree;
  theTree = 0;

  delete aFact;
  aFact = 0;
}
/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::getInstance()
{
  if (instance == 0) instance = new HadrontherapyAnalysisManager;
  return instance;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::book() 
{
  // Build up  the  analysis factory
  aFact = AIDA_createAnalysisFactory();
  AIDA::ITreeFactory* treeFact = aFact -> createTreeFactory();

  // Create the .hbk or the .root file
  G4String fileName = "DoseDistribution.hbk";
  G4String rootFileName = "DoseDistribution.root";
  
  std::string opts = "export=root";
 
  theTree = treeFact -> create(fileName,"hbook",false,true);
  theTree = treeFact -> create(rootFileName,"ROOT",false,true,opts);

  // Factories are not "managed" by an AIDA analysis system.
  // They must be deleted by the AIDA user code.
  delete treeFact;

  // Create the histogram and the ntuple factory
  histFact = aFact -> createHistogramFactory(*theTree);
  tupFact = aFact -> createTupleFactory(*theTree);

  // Create the histograms with the enrgy deposit along the X axis
  h1 = histFact -> createHistogram1D("10","slice, energy", 400, 0., 400. );

  h2 = histFact -> createHistogram1D("20","Secondary protons - slice, energy", 400, 0., 400. );
 
  h3 = histFact -> createHistogram1D("30","Secondary neutrons - slice, energy", 400, 0., 400. );

  h4 = histFact -> createHistogram1D("40","Secondary alpha - slice, energy", 400, 0., 400. );

  h5 = histFact -> createHistogram1D("50","Secondary gamma - slice, energy", 400, 0., 400. );

  h6 = histFact -> createHistogram1D("60","Secondary electron - slice, energy", 400, 0., 400. );

  h7 = histFact -> createHistogram1D("70","Secondary triton - slice, energy", 400, 0., 400. );

  h8 = histFact -> createHistogram1D("80","Secondary deuteron - slice, energy", 400, 0., 400. );

  h9 = histFact -> createHistogram1D("90","Secondary pion - slice, energy", 400, 0., 400. );
 
  h10 = histFact -> createHistogram1D("100","Energy distribution of secondary electrons", 70, 0., 70. );
 
  h11 = histFact -> createHistogram1D("110","Energy distribution of secondary photons", 70, 0., 70. );

  h12 = histFact -> createHistogram1D("120","Energy distribution of secondary deuterons", 70, 0., 70. );
 
  h13 = histFact -> createHistogram1D("130","Energy distribution of secondary tritons", 70, 0., 70. );

  h14 = histFact -> createHistogram1D("140","Energy distribution of secondary alpha particles", 70, 0., 70. );

  // Create the ntuple
  G4String columnNames = "int i; int j; int k; double energy;";
  G4String options = "";
  if (tupFact) ntuple = tupFact -> create("1","1",columnNames, options);

  // Create the ntuple
  G4String columnNames2 = "int a; double z;  int occupancy; double energy;";
  G4String options2 = "";
  if (tupFact) ionTuple = tupFact -> create("2","2", columnNames2, options2);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::FillEnergyDeposit(G4int i, 
						     G4int j, 
						     G4int k,
						     G4double energy)
{
 if (ntuple)     {
       G4int iSlice = ntuple -> findColumn("i");
       G4int jSlice = ntuple -> findColumn("j");
       G4int kSlice = ntuple -> findColumn("k");
       G4int iEnergy = ntuple -> findColumn("energy");
      
       ntuple -> fill(iSlice,i);
       ntuple -> fill(jSlice,j); 
       ntuple -> fill(kSlice,k);
       ntuple -> fill(iEnergy, energy);     }

  ntuple -> addRow(); 
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::BraggPeak(G4int slice, G4double energy)
{
  h1 -> fill(slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryProtonEnergyDeposit(G4int slice, G4double energy)
{
  h2 -> fill(slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryNeutronEnergyDeposit(G4int slice, G4double energy)
{
  h3 -> fill(slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryAlphaEnergyDeposit(G4int slice, G4double energy)
{
  h4 -> fill(slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryGammaEnergyDeposit(G4int slice, G4double energy)
{
  h5 -> fill(slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryElectronEnergyDeposit(G4int slice, G4double energy)
{
  h6 -> fill(slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryTritonEnergyDeposit(G4int slice, G4double energy)
{
  h7 -> fill(slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryDeuteronEnergyDeposit(G4int slice, G4double energy)
{
  h8 -> fill(slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::SecondaryPionEnergyDeposit(G4int slice, G4double energy)
{
  h9 -> fill(slice,energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::electronEnergyDistribution(G4double energy)
{
  h10 -> fill(energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::gammaEnergyDistribution(G4double energy)
{
  h11 -> fill(energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::deuteronEnergyDistribution(G4double energy)
{
  h12 -> fill(energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::tritonEnergyDistribution(G4double energy)
{
  h13 -> fill(energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::alphaEnergyDistribution(G4double energy)
{
  h14 -> fill(energy);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::genericIonInformation(G4int a, 
							 G4double z, 
							 G4int electronOccupancy,
							 G4double energy) 
{
  if (ionTuple)    {
       G4int aIndex = ionTuple -> findColumn("a");
       G4int zIndex = ionTuple -> findColumn("z");
       G4int electronIndex = ionTuple -> findColumn("occupancy");  
       G4int energyIndex = ionTuple -> findColumn("energy");
      
       ionTuple -> fill(aIndex,a);
      ionTuple -> fill(zIndex,z);  
      ionTuple -> fill(electronIndex, electronOccupancy); 
       ionTuple -> fill(energyIndex, energy);
     }
   ionTuple -> addRow(); 
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisManager::finish() 
{  
  // Write all histograms to file
  theTree -> commit();
 
  // Close (will again commit)
  theTree ->close();
}
#endif











