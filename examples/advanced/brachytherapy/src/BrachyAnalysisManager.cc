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
/*
Author: Susanna Guatelli
*/
// The class BrachyAnalysisManager creates and manages histograms and ntuples

// The analysis was included in this application following the extended Geant4
// example analysis/AnaEx01

#include "BrachyAnalysisManager.hh"

BrachyAnalysisManager* BrachyAnalysisManager::instance = 0;

BrachyAnalysisManager::BrachyAnalysisManager()
{
#ifdef ANALYSIS_USE
 theTFile = 0;
 histo2 =0;
 histo = 0;
#endif
}

BrachyAnalysisManager::~BrachyAnalysisManager() 
{ 
#ifdef G4ANALYSIS_USE
 delete theTFile; theTFile = 0;
 delete histo2; histo2 = 0;
 delete histo; histo = 0;
#endif
}

BrachyAnalysisManager* BrachyAnalysisManager::GetInstance()
{
  if (instance == 0) instance = new BrachyAnalysisManager;
  return instance;
}

void BrachyAnalysisManager::book() 
{  
#ifdef ANALYSIS_USE
 delete theTFile;
 theTFile = new TFile("brachytherapy.root", "RECREATE");
 
 histo = new TH1F("h10","energy spectrum", 800, 0., 800);
 histo2 = new TH2F("h20","edep2Dxy", 801, -100.125, 100.125, // binning, xmin, xmax, along x direction in mm
				     801, -100.125, 100.125);// binning, ymin, ymax, along y direction in mm
#endif 
 }

#ifdef ANALYSIS_USE
void BrachyAnalysisManager::FillH2WithEnergyDeposition(G4double xx,
                                                     G4double yy, 
                                                     G4double energyDep)
{
  histo2 -> Fill(xx, yy,energyDep);
}

void BrachyAnalysisManager::FillPrimaryParticleHistogram(G4double primaryParticleEnergy)
{
 // 1DHistogram: energy spectrum of primary particles  
  histo-> Fill(primaryParticleEnergy);
}
#endif 

void BrachyAnalysisManager::save() 
{  
#ifdef ANALYSIS_USE
 if (theTFile)
    {
	theTFile -> Write(); 
	theTFile -> Close();
    }
#endif
}
