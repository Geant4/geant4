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
#include "UltraRunActionMessenger.hh"
#include "Randomize.hh"
#include "UltraAnalysisManager.hh"

#include "Randomize.hh"


// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraRunAction::UltraRunAction()
{

   G4int seed = 1;      // RANLUX seed
   G4int luxury = 3;    // RANLUX luxury level (3 is default)
   CLHEP::HepRandom::setTheSeed(seed,luxury);
   CLHEP::HepRandom::showEngineStatus();

  theRunActMessenger = new UltraRunActionMessenger(this);
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraRunAction::~UltraRunAction()
{ 
  delete theRunActMessenger;
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraRunAction::BeginOfRunAction(const G4Run* aRun)
{

  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Open an output file
  man->OpenFile("ultra");
  man->SetFirstHistoId(1);

  // Create histogram(s)
  man->CreateH1("1","Optical photons energy (eV)", //histoID,histo name 
		500,0.,5.); //bins' number, xmin, xmax
  man->CreateH1("2","Number of Detected Photons", 
		10,0.,10.); //bins' number, xmin, xmax

  G4cout << "### Run " << aRun->GetRunID() << " start..." << G4endl;

  // save Rndm status
  if (saveRndm > 0)
    { 
      CLHEP::HepRandom::showEngineStatus();
      CLHEP::HepRandom::saveEngineStatus("beginOfRun.rndm");
    }

}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraRunAction::EndOfRunAction(const G4Run* aRun)
{

  // save Rndm status
  if (saveRndm == 1)
    { CLHEP::HepRandom::showEngineStatus();
      CLHEP::HepRandom::saveEngineStatus("endOfRun.rndm");
    }

  // Write histograms to file
  G4cout << "Close and Save Histograms" << G4endl;
  G4cout << "### Run " << aRun->GetRunID() << " ended." << G4endl;

  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  // Complete clean-up
  delete G4AnalysisManager::Instance();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraRunAction::MySetRunID(G4int myrun)
{
     runID = myrun ;
}
