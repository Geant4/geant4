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
/// \file PrimaryKiller.cc
/// \brief Implementation of the PrimaryKiller class

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

#include "PrimaryKiller.hh"

#include <G4Event.hh>
#include <G4RunManager.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UnitsTable.hh>

/** \file PrimaryKiller.cc
    \class PrimaryKiller

    Kill the primary particle:
    - either after a given energy loss
    - or after the primary particle has reached a given energy
 */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PrimaryKiller::PrimaryKiller(const G4String &name, const G4int depth)
  : G4VPrimitiveScorer(name, depth), G4UImessenger() {
  fpELossUI = std::make_unique<G4UIcmdWithADoubleAndUnit>
      ("/primaryKiller/eLossMin", this);
  fpAbortEventIfELossUpperThan = std::make_unique<G4UIcmdWithADoubleAndUnit>
      ("/primaryKiller/eLossMax", this);
  fpMinKineticE = std::make_unique<G4UIcmdWithADoubleAndUnit>
      ("/primaryKiller/minKineticE", this);
  fpSizeUI = std::make_unique<G4UIcmdWith3VectorAndUnit>
      ("/primaryKiller/setSize", this);
  fpSizeUI->SetDefaultUnit("um");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PrimaryKiller::SetNewValue(G4UIcommand *command, G4String newValue) {
  if (command == fpELossUI.get()) {
    fELossRange_Min = fpELossUI->GetNewDoubleValue(newValue);
  } else if (command == fpAbortEventIfELossUpperThan.get()) {
    fELossRange_Max = fpAbortEventIfELossUpperThan->GetNewDoubleValue(newValue);
  } else if (command == fpSizeUI.get()) {
    fPhantomSize = fpSizeUI->GetNew3VectorValue(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool PrimaryKiller::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
  const G4Track *track = aStep->GetTrack();
  G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition();

  if (std::abs(pos.x()) > fPhantomSize.getX() / 2 || std::abs(pos.y()) > fPhantomSize.getY() / 2
      || std::abs(pos.z()) > fPhantomSize.getZ() / 2) {
    const_cast<G4Track *>(track)->SetTrackStatus(fStopAndKill);
    return false;
  }

  if (track->GetTrackID() != 1 || track->GetParticleDefinition()->GetPDGEncoding() != 11) {
    return FALSE;
  }

  //-------------------

  const G4double kineticE = aStep->GetPostStepPoint()->GetKineticEnergy();

  const G4double eLoss = aStep->GetPreStepPoint()->GetKineticEnergy() - kineticE;

  if (eLoss == 0.) { return FALSE; }

  //-------------------

  fELoss += eLoss;

  if (fELoss > fELossRange_Max) {
    G4RunManager::GetRunManager()->AbortEvent();
    /*    int eventID =
         G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

        G4cout << " * PrimaryKiller: aborts event " << eventID <<" energy loss "
                  "is too large. \n"
               << " * Energy loss by primary is: "
               << G4BestUnit(fELoss, "Energy")
               << ". Event is aborted if the Eloss is > "
               << G4BestUnit(fELossRange_Max, "Energy")
               << G4endl;
     */
  }

  if (fELoss >= fELossRange_Min || kineticE <= fKineticE_Min) {
    const_cast<G4Track *>(track)->SetTrackStatus(fStopAndKill);
    //     G4cout << "kill track at : "<<'\n';
    //           << G4BestUnit(kineticE, "Energy")
    //           << ", E loss is: "
    //           << G4BestUnit(fELoss, "Energy")
    //           << " /fELossMax: "
    //           << G4BestUnit(fELossMax, "Energy")
    //           << ", EThreshold is: "
    //           << G4BestUnit(fEThreshold, "Energy")
    //           << G4endl;
  }

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
