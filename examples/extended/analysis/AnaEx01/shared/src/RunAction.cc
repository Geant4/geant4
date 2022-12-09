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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4AccumulableManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(HistoManager* histo)
: fHistoManager(histo)
{
  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fSumEAbs);
  accumulableManager->RegisterAccumulable(fSum2EAbs);
  accumulableManager->RegisterAccumulable(fSumEGap);
  accumulableManager->RegisterAccumulable(fSum2EGap);
  accumulableManager->RegisterAccumulable(fSumLAbs);
  accumulableManager->RegisterAccumulable(fSum2LAbs);
  accumulableManager->RegisterAccumulable(fSumLGap);
  accumulableManager->RegisterAccumulable(fSum2LGap);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // reset accumulables to their initial values
  G4AccumulableManager::Instance()->Reset();

  //histograms
  //
  fHistoManager->Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillPerEvent(G4double EAbs, G4double EGap,
                             G4double LAbs, G4double LGap)
{
  //accumulate statistic
  //
  fSumEAbs += EAbs;  fSum2EAbs += EAbs*EAbs;
  fSumEGap += EGap;  fSum2EGap += EGap*EGap;

  fSumLAbs += LAbs;  fSum2LAbs += LAbs*LAbs;
  fSumLGap += LGap;  fSum2LGap += LGap*LGap;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  // Merge accumulables
  G4AccumulableManager::Instance()->Merge();

  G4int nofEvents = aRun->GetNumberOfEvent();
  if (nofEvents == 0) {
    // close open files
    fHistoManager->Save();
    return;
  }

  //compute statistics: mean and rms
  //
  auto sumEAbs = fSumEAbs.GetValue();
  auto sum2EAbs = fSum2EAbs.GetValue();
  sumEAbs /= nofEvents; sum2EAbs /= nofEvents;
  auto rmsEAbs = sum2EAbs - sumEAbs*sumEAbs;
  if (rmsEAbs >0.) {
    rmsEAbs = std::sqrt(rmsEAbs);
  } else {
    rmsEAbs = 0.;
  }

  auto sumEGap = fSumEGap.GetValue();
  auto sum2EGap = fSum2EGap.GetValue();
  sumEGap /= nofEvents; sum2EGap /= nofEvents;
  auto rmsEGap = sum2EGap - sumEGap*sumEGap;
  if (rmsEGap >0.) {
    rmsEGap = std::sqrt(rmsEGap);
  } else {
    rmsEGap = 0.;
  }

  auto sumLAbs = fSumLAbs.GetValue();
  auto sum2LAbs = fSum2LAbs.GetValue();
  sumLAbs /= nofEvents; sum2LAbs /= nofEvents;
  auto rmsLAbs = sum2LAbs - sumLAbs*sumLAbs;
  if (rmsLAbs >0.) {
    rmsLAbs = std::sqrt(rmsLAbs);
  } else {
    rmsLAbs = 0.;
  }

  auto sumLGap = fSumLGap.GetValue();
  auto sum2LGap = fSum2LGap.GetValue();
  sumLGap /= nofEvents; sum2LGap /= nofEvents;
  G4double rmsLGap = sum2LGap - sumLGap*sumLGap;
  if (rmsLGap >0.) {
    rmsLGap = std::sqrt(rmsLGap);
  } else {
    rmsLGap = 0.;
  }

  //print
  //
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in Absorber : " << G4BestUnit(sumEAbs,"Energy")
     << " +- "                          << G4BestUnit(rmsEAbs,"Energy")
     << "\n mean Energy in Gap      : " << G4BestUnit(sumEGap,"Energy")
     << " +- "                          << G4BestUnit(rmsEGap,"Energy")
     << G4endl;

  G4cout
     << "\n mean trackLength in Absorber : " << G4BestUnit(sumLAbs,"Length")
     << " +- "                               << G4BestUnit(rmsLAbs,"Length")
     << "\n mean trackLength in Gap      : " << G4BestUnit(sumLGap,"Length")
     << " +- "                               << G4BestUnit(rmsLGap,"Length")
     << "\n------------------------------------------------------------\n"
     << G4endl;

  //save histograms
  //
  fHistoManager->PrintStatistic();
  fHistoManager->Save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
