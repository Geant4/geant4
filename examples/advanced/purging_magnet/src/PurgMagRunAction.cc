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
//    ******************************
//    *                            *
//    *    PurgMagRunAction.cc     *
//    *                            *
//    ******************************
//
//

#include "PurgMagRunAction.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include "PurgMagAnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagRunAction::PurgMagRunAction()
{   
  saveRndm = 1;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagRunAction::~PurgMagRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagRunAction::BeginOfRunAction(const G4Run* aRun)
{  
  //Analysis must be handled by master and workers
  Book();  
  
  if (IsMaster())    
    G4cout << "---> Run " << aRun->GetRunID() << " (master) start." 
	   << G4endl;
  else
    G4cout << "---> Run " << aRun->GetRunID() << " (worker) start." 
	   << G4endl;

  
  // save Rndm status
  if (saveRndm > 0 && IsMaster())
    { 
      G4Random::showEngineStatus();
      G4Random::saveEngineStatus("beginOfRun.rndm");
    }
      
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagRunAction::EndOfRunAction(const G4Run* aRun)
{      
  if (IsMaster())    
    G4cout << "Total number of event = " << aRun->GetNumberOfEvent() << G4endl;
  else
    G4cout << "Partial number of event in this worker = " 
	   << aRun->GetNumberOfEvent() << G4endl;
 
       
  if (IsMaster())
    {
      // save Rndm status
      if (saveRndm == 1)
	{ 
	  G4Random::showEngineStatus();
	  G4Random::saveEngineStatus("endOfRun.rndm");
	}             
    }            
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagRunAction::Book() 
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



