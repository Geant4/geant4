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

<<<<<<< HEAD
#include <stdlib.h>
#include "BrachyAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
=======
#include "BrachyAnalysisManager.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

BrachyAnalysisManager* BrachyAnalysisManager::instance = 0;

BrachyAnalysisManager::BrachyAnalysisManager()
{
#ifdef ANALYSIS_USE
 theTFile = 0;
 histo = 0;
 ntuple = 0; 
#endif
}

BrachyAnalysisManager::~BrachyAnalysisManager() 
{ 
#ifdef G4ANALYSIS_USE
 delete theTFile; theTFile = 0;
 delete histo; histo = 0;
 delete ntuple; ntuple = 0; 
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
 ntuple =  new TNtuple("ntuple","edep3D","xx:yy:zz:edep");
#endif 
 }

#ifdef ANALYSIS_USE
void BrachyAnalysisManager::FillNtupleWithEnergyDeposition(G4double xx,
                                                     G4double yy, 
                                                     G4double zz,
                                                     G4double energyDep)
{
  ntuple -> Fill(xx, yy, zz, energyDep);
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
