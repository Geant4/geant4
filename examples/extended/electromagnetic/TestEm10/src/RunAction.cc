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
/// \file electromagnetic/TestEm10/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
//
// 

#include "RunAction.hh"
#include "RunMessenger.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4ios.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction()
  : G4UserRunAction(),
    fRunMessenger(0), fRndmFreq(1)
{
  fRunMessenger = new RunMessenger(this);

  BookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BookHisto()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName("testem10");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);     // start histogram numbering from 1
  analysisManager->SetActivation(true);    // enable inactivation of histograms
  
  // Define histograms start values
  const G4int kMaxHisto = 5;
  const G4String id[] = 
                { "1",  "2", "3", "4", "5"};
  const G4String title[] = 
                { "Edep",                      // 1
                  "XTR Gamma spectrum",        // 2
                  "Secondary Gamma spectrum" , // 3
                  "Secondary e- spectrum",     // 4
                  "Edep.old"};                 // 5
   // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  if (fRndmFreq > 0) { 
    CLHEP::HepRandom::showEngineStatus();
    CLHEP::HepRandom::saveEngineStatus("beginOfRun.rndm");
  }  

  //  histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run* run)
{
  // print run statisctics
  G4int nofEvents = run->GetNumberOfEvent();
  G4cout << " ================== run summary =====================" << G4endl;
  G4cout << " End of Run TotNbofEvents = "  << nofEvents << G4endl ;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << " Mean energy deposit in absorber = "
           << analysisManager->GetH1(1)->mean()/MeV << " +-"
           << analysisManager->GetH1(1)->rms()/MeV  << "  MeV " << G4endl;
  }
  if ( analysisManager->GetH1(2) ) {
    G4cout << " Total number of XTR gammas = "
           << analysisManager->GetH1(2)->entries() << G4endl;
  }
  if ( analysisManager->GetH1(3) ) {
    G4cout << " Total number of all gammas = "
           << analysisManager->GetH1(3)->entries() << G4endl;
  }

  // save Rndm status
  if (fRndmFreq == 1) { 
    CLHEP::HepRandom::showEngineStatus();
    CLHEP::HepRandom::saveEngineStatus("endOfRun.rndm");
  }

  //save histograms      
  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
