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
/// \file analysis/N03Con/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ConvergenceTester.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    
  //initialize cumulative quantities
  //
  fSumEAbs = fSum2EAbs = fSumEGap = fSum2EGap = 0.;
  fSumLAbs = fSum2LAbs = fSumLGap = fSum2LGap = 0.;

  fEabs_tally = new G4ConvergenceTester();
  fEgap_tally = new G4ConvergenceTester();
  fLabs_tally = new G4ConvergenceTester();
  fLgap_tally = new G4ConvergenceTester();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillPerEvent(G4double eAbs, G4double eGap,
                                  G4double lAbs, G4double lGap)
{
  //accumulate statistic
  //
  fSumEAbs += eAbs;  fSum2EAbs += eAbs*eAbs;
  fSumEGap += eGap;  fSum2EGap += eGap*eGap;
  
  fSumLAbs += lAbs;  fSum2LAbs += lAbs*lAbs;
  fSumLGap += lGap;  fSum2LGap += lGap*lGap;

  fEabs_tally->AddScore( eAbs ); 
  fEgap_tally->AddScore( eGap ); 
  fLabs_tally->AddScore( lAbs ); 
  fLgap_tally->AddScore( lGap );     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  //compute statistics: mean and rms
  //
  fSumEAbs /= NbOfEvents; fSum2EAbs /= NbOfEvents;
  G4double rmsEAbs = fSum2EAbs - fSumEAbs*fSumEAbs;
  if (rmsEAbs >0.) rmsEAbs = std::sqrt(rmsEAbs); else rmsEAbs = 0.;
  
  fSumEGap /= NbOfEvents; fSum2EGap /= NbOfEvents;
  G4double rmsEGap = fSum2EGap - fSumEGap*fSumEGap;
  if (rmsEGap >0.) rmsEGap = std::sqrt(rmsEGap); else rmsEGap = 0.;
  
  fSumLAbs /= NbOfEvents; fSum2LAbs /= NbOfEvents;
  G4double rmsLAbs = fSum2LAbs - fSumLAbs*fSumLAbs;
  if (rmsLAbs >0.) rmsLAbs = std::sqrt(rmsLAbs); else rmsLAbs = 0.;
  
  fSumLGap /= NbOfEvents; fSum2LGap /= NbOfEvents;
  G4double rmsLGap = fSum2LGap - fSumLGap*fSumLGap;
  if (rmsLGap >0.) rmsLGap = std::sqrt(rmsLGap); else rmsLGap = 0.;
  
  //print
  //
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in Absorber : " << G4BestUnit(fSumEAbs,"Energy")
     << " +- "                          << G4BestUnit(rmsEAbs,"Energy")  
     << "\n mean Energy in Gap      : " << G4BestUnit(fSumEGap,"Energy")
     << " +- "                          << G4BestUnit(rmsEGap,"Energy")
     << G4endl;
     
  G4cout
     << "\n mean trackLength in Absorber : " << G4BestUnit(fSumLAbs,"Length")
     << " +- "                               << G4BestUnit(rmsLAbs,"Length")  
     << "\n mean trackLength in Gap      : " << G4BestUnit(fSumLGap,"Length")
     << " +- "                               << G4BestUnit(rmsLGap,"Length")
     << "\n------------------------------------------------------------\n"
     << G4endl;

  fEabs_tally->ShowResult();
  fEabs_tally->ShowHistory();
  fEgap_tally->ShowResult();
  fEgap_tally->ShowHistory();

  fLabs_tally->ShowResult();
  fLabs_tally->ShowHistory();
  fLgap_tally->ShowResult();
  fLgap_tally->ShowHistory();     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
