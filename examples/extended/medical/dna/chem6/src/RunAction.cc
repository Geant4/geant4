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
// chem6 example is derived from chem4 and chem5 examples
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
// Authors: W. G. Shin and S. Incerti (CENBG, France)
//
// $Id$
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

extern std::ofstream out;

RunAction::RunAction()
 : G4UserRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4Run* RunAction::GenerateRun()
{
  Run* run = new Run();
  return run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " starts." << G4endl;

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
  const Run* chem6Run = static_cast<const Run*>(run);
  G4double sumDose   = chem6Run->GetSumDose();

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
      dynamic_cast<ScoreSpecies*>(chem6Run->GetPrimitiveScorer());

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

    if(nOfEvent > 1)
    {
      out<<std::setw(12)<<"LET"<<std::setw(12)<<LET_mean
          <<std::setw(12)<<"LET_SD"<<std::setw(12)<<LET_square/(nOfEvent-1)<<'\n';
    }else
    {
      out<<std::setw(12)<<"LET"<<std::setw(12)<<LET_mean
          <<std::setw(12)<<"LET_SD"<<std::setw(12)<<LET_square<<'\n';
    }

    masterScorer->OutputAndClear();

    out<<'\n';

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
    << " Total energy deposited in the world volume : " << sumDose/eV << " eV"
    << G4endl
    << " ------------------------------------------------------------"
    << G4endl
    << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
