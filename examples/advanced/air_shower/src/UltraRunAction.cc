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
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraRunAction.cc
//    ****************************************************
//
//    RunAction class for Ultra; it has also a Messenger Class
//
#include "UltraRunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "UltraAnalysisManager.hh"
#include <ctime>
#include "Randomize.hh"


// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraRunAction::UltraRunAction()
{
   seed = -1;      // RANLUX seed
   luxury = 3;     // RANLUX luxury level (3 is default)
   saveRndm = 1;
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraRunAction::~UltraRunAction()
{;}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraRunAction::BeginOfRunAction(const G4Run* aRun)
{
  // Get/create analysis manager: need to do that in the master and in the workers
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Open an output file
  man->OpenFile("ultra");
  man->SetFirstHistoId(1);

  // Create histogram(s)
  man->CreateH1("1","Optical photons energy (eV)", //histoID,histo name 
		500,0.,5.); //bins' number, xmin, xmax
  man->CreateH1("2","Number of Detected Photons", 
		10,0.,10.); //bins' number, xmin, xmax

  if (!IsMaster()) //it is a slave, do nothing else
    {
      G4cout << "ooo Run " << aRun->GetRunID() << " starts on slave." << G4endl;
      return;
    }

  //Master or sequential
  G4cout << "ooo Run " << aRun->GetRunID() << " starts (global)." << G4endl;
  if (seed<0) //not initialized by anybody else
    {
      seed=time(0);
      G4Random::setTheSeed(seed,luxury);
      G4Random::showEngineStatus();
    }      
      
  // save Rndm status
  if (saveRndm > 0)       
    G4Random::saveEngineStatus("beginOfRun.rndm");

  return; 
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraRunAction::EndOfRunAction(const G4Run* aRun)
{
  // Write histograms to file  
  

  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  // Complete clean-up
  delete G4AnalysisManager::Instance();

  if (!IsMaster())
    {
      G4cout << "### Run " << aRun->GetRunID() << " (slave) ended." << G4endl;
      return;
    }

  G4cout << "### Run " << aRun->GetRunID() << " (global) ended." << G4endl;
  // save Rndm status
  if (saveRndm == 1)        
    G4Random::saveEngineStatus("endOfRun.rndm");
  return;
}


