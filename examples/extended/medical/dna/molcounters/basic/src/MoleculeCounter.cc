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
/// \file MoleculeCounter.cc
/// \brief Implementation of the MoleculeCounter class

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025)
//

#include <memory>

#include "MoleculeCounter.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4MoleculeLocator.hh"
#include "G4Navigator.hh"
#include "G4Scheduler.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
using G4VMoleculeCounterIndex = G4VMoleculeCounter::G4VMoleculeCounterIndex;

G4ThreadLocal std::unique_ptr<G4Navigator> MoleculeCounter::fNavigator = nullptr;

MoleculeCounter::MoleculeCounter(G4String name)
  : G4VUserMoleculeCounter(name, G4VMoleculeCounter::MoleculeCounterType::Other),
    fIgnoreMoleculePosition(false),
    fIgnoreCopyNumbers(true)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MoleculeCounter::InitializeUser() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::unique_ptr<G4VMoleculeCounterIndex> MoleculeCounter::BuildIndex(const G4Track* aTrack) const
{
  if (fVerbose > 1) {
    G4cout << "MoleculeCounter::BuildIndex(" << aTrack->GetTrackID() << " : "
           << GetMolecule(aTrack)->GetName() << ")" << G4endl;
  }
  if (fIgnoreMoleculePosition) {
    return std::make_unique<MoleculeCounterIndex>(
      GetMolecule(aTrack)->GetMolecularConfiguration(), nullptr);
  }
  else {
    const G4VTouchable* touchable = aTrack->GetNextTouchable();
    G4TouchableHistory* touchableHistory = nullptr;

    if (touchable == nullptr) touchable = aTrack->GetTouchable();
    if (touchable == nullptr) {
      auto touchableHandle = G4MoleculeLocator::Instance()->LocateMoleculeTrack(aTrack);
      touchable = touchableHandle();
    }
    if (touchable == nullptr) {  // still not found -- should never fire
      G4ExceptionDescription errMsg;
      errMsg << "Molecule scorer requires a valid volume pointer."
             << " G4Track->GetVolume: " << aTrack->GetVolume()
             << " G4Navigator->LocateGlobalPointAndUpdateTouchable: nullptr" << G4endl;
      G4Exception("MoleculeCounter::BuildIndex", "VOL_NOT_FOUND", FatalException, errMsg);
    }

    auto index = std::make_unique<MoleculeCounterIndex>(
        GetMolecule(aTrack)->GetMolecularConfiguration(), touchable);

    delete touchableHistory;
    return index;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::unique_ptr<G4VMoleculeCounterIndex>
MoleculeCounter::BuildIndex(const G4Track* aTrack, const G4StepPoint* aStepPoint) const
{
  return std::make_unique<MoleculeCounterIndex>(
      GetMolecule(aTrack)->GetMolecularConfiguration(), aStepPoint->GetTouchable());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::unique_ptr<G4VMoleculeCounterIndex>
MoleculeCounter::BuildSimpleIndex(const G4MolecularConfiguration* configuration) const
{
  return std::make_unique<MoleculeCounterIndex>(
    configuration, nullptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool MoleculeCounter::GetIgnoreMoleculePosition() const
{
  return fIgnoreMoleculePosition;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MoleculeCounter::SetIgnoreMoleculePosition(G4bool flag)
{
  fIgnoreMoleculePosition = flag;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MoleculeCounter::SetNegativeCountsAreFatal(G4bool flag)
{
  G4VMoleculeCounter::SetNegativeCountsAreFatal(flag);
  // this is protected in G4VMoleculeCounter and thus, must be exposed by the derived class
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
