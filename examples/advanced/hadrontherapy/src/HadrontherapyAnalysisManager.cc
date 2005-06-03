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
// $Id: HadrontherapyAnalisysManager.cc;  May 2005
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

HadrontherapyAnalysisManager::HadrontherapyAnalysisManager() : 
  aFact(0), theTree(0), histFact(0), tupFact(0), h1(0), ntuple(0)
{  
}

HadrontherapyAnalysisManager::~HadrontherapyAnalysisManager() 
{ 
  delete ntuple;
  ntuple = 0;

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

HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::getInstance()
{
  if (instance == 0) instance = new HadrontherapyAnalysisManager;
  return instance;
}

void HadrontherapyAnalysisManager::book() 
{
  // Build up  the  analysis factory
  aFact = AIDA_createAnalysisFactory();
  AIDA::ITreeFactory* treeFact = aFact -> createTreeFactory();

  // Create the .hbk file
  G4String fileName = "hadrontherapy.hbk";
  theTree = treeFact -> create(fileName,"hbook",false,true);
  delete treeFact;

  // Create the histogram and the ntuple factory
  histFact = aFact -> createHistogramFactory(*theTree);
  tupFact = aFact -> createTupleFactory(*theTree);

  // Create the histogram
  h1 = histFact -> createHistogram1D("10","slice, energy", 80, 0., 80. );

  // Create the ntuple
  G4String columnNames = "int i; int j; int k; double energy;";
  G4String options = "";
  if (tupFact) ntuple = tupFact->create("1","1",columnNames, options);
}

void HadrontherapyAnalysisManager::FillEnergyDeposit(G4int i, 
						     G4int j, 
						     G4int k,
						     G4double energy)
{
  if (ntuple)
    {
      G4int iSlice = ntuple -> findColumn("i");
      G4int jSlice = ntuple -> findColumn("j");
      G4int kSlice = ntuple -> findColumn("k");
      G4int iEnergy = ntuple -> findColumn("energy");
      
      ntuple -> fill(iSlice,i);
      ntuple -> fill(jSlice,j); 
      ntuple -> fill(kSlice,k);
      ntuple -> fill(iEnergy, energy);
    }

  ntuple -> addRow(); 
}

void HadrontherapyAnalysisManager::BraggPeak(G4int slice, G4double energy)
{
  h1 -> fill(slice,energy);
}

void HadrontherapyAnalysisManager::finish() 
{  
 // Write all histograms to file
 theTree -> commit();
 
 // Close (will again commit)
 theTree ->close();
}
#endif











