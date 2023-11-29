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
// This example is provided by the Geant4-DNA collaboration
// dnadamage3 example is derived from the chem6 example
// chem6 example authors: W. G. Shin and S. Incerti (CENBG, France)
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Run.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

RunAction::RunAction()
 : G4UserRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4Run* RunAction::GenerateRun()
{
  Run* run = new Run();
  return run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::BeginOfRunAction(const G4Run*)
{
  // informs the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // results
  //
  const Run* dnadamage3Run = static_cast<const Run*>(run);
  G4double sumDose   = dnadamage3Run->GetSumDose();

  // print
  //
  if (IsMaster())
  {

    G4cout
      << G4endl
      << "--------------------End of Global Run-----------------------"
      << G4endl
      << "  The run has " << nofEvents << " events "
      << G4endl;

    ScoreSpecies* masterScorer=
      dynamic_cast<ScoreSpecies*>(dnadamage3Run->GetPrimitiveScorer());

    ScoreStrandBreaks* masterSBScorer=
      dynamic_cast<ScoreStrandBreaks*>(dnadamage3Run->GetSBScorer());

    G4cout << "Number of events recorded by the species scorer = "
      << masterScorer->GetNumberOfRecordedEvents()
      << G4endl;

    // LET
    Run* aRun = (Run*)run;
    G4THitsMap<G4double>* totLET = aRun->GetLET();
    G4int nOfEvent = totLET->entries();
    G4double LET_mean = 0;
    G4double LET_square = 0;
    for(G4int i=0;i<nOfEvent;i++){
      G4double* LET = (*totLET)[i];
      if(!LET) continue;
      LET_mean += *LET;
      LET_square += (*LET)*(*LET);
    }
    LET_mean /= nOfEvent;
    LET_square = std::sqrt(LET_square/nOfEvent - std::pow(LET_mean,2));

    masterScorer->OutputAndClear();
    masterSBScorer->OutputAndClear(LET_mean,LET_square);

  }
  else
  {
    G4cout
      << G4endl
      << "--------------------End of Local Run------------------------"
      << G4endl
      << "  The run has " << nofEvents << " events"
      << G4endl;
  }

  G4cout
    << " Total energy deposited in the world volume : "
    << sumDose/eV << " eV"
    << G4endl
    << " ------------------------------------------------------------"
    << G4endl
    << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
