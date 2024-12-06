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
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file TimeStepAction.hh
/// \brief Implementation of the TimeStepAction class

#include "TimeStepAction.hh"

#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4MoleculeTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <G4Scheduler.hh>

TimeStepAction::TimeStepAction() : G4UserTimeStepAction()
{
  /**
   * Give to G4ITTimeStepper the user defined time steps
   * eg : from 1 picosecond to 10 picosecond, the minimum time
   * step that the TimeStepper can returned is 0.1 picosecond.
   * Those time steps are used for the chemistry of G4DNA
   */

  AddTimeStep(1 * picosecond, 0.1 * picosecond);
  AddTimeStep(10 * picosecond, 1 * picosecond);
  AddTimeStep(100 * picosecond, 3 * picosecond);
  AddTimeStep(1000 * picosecond, 10 * picosecond);
  AddTimeStep(10000 * picosecond, 100 * picosecond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserPostTimeStepAction()
{
  if (G4Scheduler::Instance()->GetGlobalTime() > 99 * ns) {
    G4cout << "_________________" << G4endl;
    G4cout << "At : " << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time") << G4endl;

    auto species = G4MoleculeTable::Instance()->GetConfiguration("°OH");
    PrintSpecieInfo(species);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserReactionAction(const G4Track& reactantA, const G4Track& reactantB,
                                        const std::vector<G4Track*>* products)
{
  // this function shows how to get species ID, positions of reaction.
  G4cout << G4endl;

  G4cout << "At : " << G4Scheduler::Instance()->GetGlobalTime() / ns
         << " (ns) reactantA = " << GetMolecule(reactantA)->GetName()
         << " (ID number = " << reactantA.GetTrackID() << ")"
         << " at position : " << reactantA.GetPosition() / nm
         << " reacts with reactantB = " << GetMolecule(&reactantB)->GetName()
         << " (ID number = " << reactantB.GetTrackID() << ")"
         << " at position : " << reactantA.GetPosition() / nm << G4endl;

  if (products) {
    auto nbProducts = (G4int)products->size();
    for (G4int i = 0; i < nbProducts; i++) {
      G4cout << "      creating product " << i + 1 << " =" << GetMolecule((*products)[i])->GetName()
             << " position : " << (*products)[i]->GetPosition() << G4endl;
    }
  }
  G4cout << G4endl;
}

void TimeStepAction::PrintSpecieInfo(G4MolecularConfiguration* molconf)
{
  // this function shows how to get a specific species ID, positions at each time step.
  auto moleculeID = molconf->GetMoleculeID();
  const G4String& moleculeName = molconf->GetFormatedName();
  G4cout << "Get inf of : " << moleculeName << G4endl;
  auto trackList = G4ITTrackHolder::Instance()->GetMainList(moleculeID);

  if (trackList == nullptr) {
    G4cout << "No species" << G4endl;
    return;
  }
  G4TrackList::iterator it = trackList->begin();
  G4TrackList::iterator end = trackList->end();
  for (; it != end; ++it) {
    auto track = *it;
    G4cout << "TrackID: " << track->GetTrackID() << "  position : " << track->GetPosition() / nm
           << G4endl;
  }
}
