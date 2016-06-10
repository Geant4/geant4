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
// Code developed by:
//  S.Larsson
//
//    ************************************
//    *                                  *
//    *    PurgMagAnalysisManager.cc     *
//    *                                  *
//    ************************************
//
// $Id: PurgMagAnalysisManager.cc 84477 2014-10-16 08:44:04Z gcosmo $
//
#include "PurgMagAnalysisManager.hh"


PurgMagAnalysisManager* PurgMagAnalysisManager::instance = 0;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagAnalysisManager::PurgMagAnalysisManager() 
{;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagAnalysisManager::~PurgMagAnalysisManager() 
{;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagAnalysisManager* PurgMagAnalysisManager::getInstance()
{
  if (instance == 0) instance = new PurgMagAnalysisManager;
  return instance;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagAnalysisManager::book() 
{
  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();
 
  // Open an output file
  man->OpenFile("purgmag");
  man->SetFirstNtupleId(1);
  
  // Create 1st ntuple (id = 1)
  //    
  man->CreateNtuple("n101", "Electron");
  man->CreateNtupleDColumn("ex");
  man->CreateNtupleDColumn("ey");
  man->CreateNtupleDColumn("ez");
  man->CreateNtupleDColumn("ee");
  man->CreateNtupleDColumn("epx");
  man->CreateNtupleDColumn("epy");
  man->CreateNtupleDColumn("epz");
  man->FinishNtuple();
  G4cout << "Ntuple-1 created" << G4endl;

  // Create 2nd ntuple (id = 2)
  //    
  man->CreateNtuple("n102", "Gamma");
  man->CreateNtupleDColumn("gx");
  man->CreateNtupleDColumn("gy");
  man->CreateNtupleDColumn("gz");
  man->CreateNtupleDColumn("ge");
  man->CreateNtupleDColumn("gpx");
  man->CreateNtupleDColumn("gpy");
  man->CreateNtupleDColumn("gpz");
  man->FinishNtuple();
  G4cout << "Ntuple-2 created" << G4endl;
 
  // Create 3rd ntuple (id = 3)
  //
  man->CreateNtuple("n103", "Positron");
  man->CreateNtupleDColumn("px");
  man->CreateNtupleDColumn("py");
  man->CreateNtupleDColumn("pz");
  man->CreateNtupleDColumn("pe");
  man->CreateNtupleDColumn("ppx");
  man->CreateNtupleDColumn("ppy");
  man->CreateNtupleDColumn("ppz");
  man->FinishNtuple();
  G4cout << "Ntuple-3 created" << G4endl;

  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fills a N-tuple with position, energy and momentum of 
// electrons entering the measurement volume. 
void PurgMagAnalysisManager::fill_Tuple_Electrons(G4double ex, G4double ey, G4double ez,     // Position
						  G4double ee,                               // Energy
						  G4double epx, G4double epy, G4double epz)  // Momentum
{
  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Fill 1st ntuple ( id = 1)
  man->FillNtupleDColumn(1,0, ex);
  man->FillNtupleDColumn(1,1, ey);
  man->FillNtupleDColumn(1,2, ez);
  man->FillNtupleDColumn(1,3, ee);
  man->FillNtupleDColumn(1,4, epx);
  man->FillNtupleDColumn(1,5, epy);
  man->FillNtupleDColumn(1,6, epz);
  man->AddNtupleRow(1);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fills a N-tuple with position, energy and momentum of 
// photons entering the measurement volume. 
void PurgMagAnalysisManager::fill_Tuple_Gamma(G4double gx, G4double gy, G4double gz,     // Position 
					      G4double ge,                               // Energy
       					      G4double gpx, G4double gpy, G4double gpz)  // Momentum
{ 
  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Fill 2nd ntuple ( id = 2)
  man->FillNtupleDColumn(2,0, gx);
  man->FillNtupleDColumn(2,1, gy);
  man->FillNtupleDColumn(2,2, gz);
  man->FillNtupleDColumn(2,3, ge);
  man->FillNtupleDColumn(2,4, gpx);
  man->FillNtupleDColumn(2,5, gpy);
  man->FillNtupleDColumn(2,6, gpz);
  man->AddNtupleRow(2);  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fills a N-tuple with position, energy and momentum of 
// positrons entering the measurement volume. 
void PurgMagAnalysisManager::fill_Tuple_Positrons(G4double px, G4double py, G4double pz,     // Position 
					     G4double pe,                               // Energy
					     G4double ppx, G4double ppy, G4double ppz)  // Momentum
{
  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Fill 3rd ntuple ( id = 3)
  man->FillNtupleDColumn(3,0, px);
  man->FillNtupleDColumn(3,1, py);
  man->FillNtupleDColumn(3,2, pz);
  man->FillNtupleDColumn(3,3, pe);
  man->FillNtupleDColumn(3,4, ppx);
  man->FillNtupleDColumn(3,5, ppy);
  man->FillNtupleDColumn(3,6, ppz);
  man->AddNtupleRow(3);  

}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void PurgMagAnalysisManager::finish() 
{  
  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  // Complete clean-up
  delete G4AnalysisManager::Instance();
}


