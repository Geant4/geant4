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
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 51 (2024) 5873-5889
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"
#include "Run.hh"

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Tubs.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4ParticleDefinition* partDef =
    step->GetTrack()->GetDynamicParticle()->GetDefinition();

  // Boolean to force kinetic energy value
  if (energyFix)
  {
    G4DNAGenericIonsManager* instance = G4DNAGenericIonsManager::Instance();

    if ( partDef == G4Proton::ProtonDefinition()
       ||
       partDef == G4Alpha::AlphaDefinition()
       ||
       partDef == instance->GetIon("hydrogen")
       ||
       partDef == instance->GetIon("alpha+")
       ||
       partDef == instance->GetIon("helium")
       ||
       partDef->GetPDGCharge()>4 )
    {
      const PrimaryGeneratorAction* primaryGenerator =
        static_cast<const PrimaryGeneratorAction*>
        (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

      G4double primaryEnergy = primaryGenerator->GetPrimaryKineticEnergy();

      // *** WARNING ***
      // this is a non physical trick to fix the kinetic energy
      // of the ion at post step
      step->GetPostStepPoint()->SetKineticEnergy(primaryEnergy);

      //G4cout << "-------------- NEW STEP ------------"
      //  <<  step->GetPreStepPoint()->GetKineticEnergy()/eV << G4endl;
      //G4cout << "-----Kinetic energy at PreStepPoint="
      //  <<  step->GetPreStepPoint()->GetKineticEnergy()/eV << G4endl;
      //G4cout << "-----Kinetic energy at PostStepPoint="
      // << step->GetPostStepPoint()->GetKineticEnergy()/eV << G4endl;
    }
  }

  if (!step->GetPostStepPoint()) return;

  if (!step->GetPostStepPoint()->GetProcessDefinedStep()) return;

  // Absorbed dose deposited by electrons only

  if (partDef == G4Electron::ElectronDefinition()) {

    G4VPhysicalVolume* volume =
      step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (!volume) return;

    G4LogicalVolume* logicVol = volume->GetLogicalVolume();
    G4double dose = step->GetTotalEnergyDeposit()/logicVol->GetMass();

    Run* run
       = static_cast<Run*>
       (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

    if (dose>0) run->AddCylDoseDeposit(volume->GetCopyNo(), dose);

    /*
    // Mass computation check

    G4VSolid* solidVol = logicVol->GetSolid();
    G4Tubs* cylinder = dynamic_cast<G4Tubs*>(solidVol);
    if (cylinder)
    {
      G4double rInner = cylinder->GetInnerRadius();
      G4double rOuter = cylinder->GetOuterRadius();
      G4cout << "**** MASS CHECK ****" << G4endl;
      G4cout << "Copy number=" << volume->GetCopyNo() << G4endl;
      G4cout << "rInner/nm=" << rInner/nm << G4endl;
      G4cout << "rOuter/nm=" << rOuter/nm << G4endl;
      G4cout << "length/nm=" << 2*cylinder->GetZHalfLength()/nm << G4endl;
      G4cout << "mass/kg=" << logicVol->GetMass()/kg << G4endl;
    }
    */
  } // electron

}
