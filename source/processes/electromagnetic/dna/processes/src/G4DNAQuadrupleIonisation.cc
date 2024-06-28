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
// G4DNAQuadrupleIonisation.cc
//
//  Created at 2024/04/03 (Thu.)
//  Author: Shogo OKADA @KEK-CRC (shogo.okada@kek.jp)
//

#include "G4DNAQuadrupleIonisation.hh"
#include "G4DNAQuadrupleIonisationModel.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericIon.hh"
#include "G4Positron.hh"

//------------------------------------------------------------------------------
G4DNAQuadrupleIonisation::G4DNAQuadrupleIonisation(
  const G4String& pname, G4ProcessType type)
    : G4VEmProcess(pname, type),
      is_initialized_(false)
{
  SetProcessSubType(fLowEnergyQuadrupleIonisation);
}

//------------------------------------------------------------------------------
G4bool G4DNAQuadrupleIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (
    &p == G4Proton::Proton() ||
    &p == G4DNAGenericIonsManager::Instance()->GetIon("alpha++") ||
    &p == G4GenericIon::GenericIonDefinition()
  );
}

//------------------------------------------------------------------------------
void G4DNAQuadrupleIonisation::InitialiseProcess(const G4ParticleDefinition* p)
{
  if (is_initialized_) { return; }

  is_initialized_ = true;
  SetBuildTableFlag(false);

  const auto name = p->GetParticleName();

  if (name == "proton") {

    if (!EmModel()) {
      auto ptr = new G4DNAQuadrupleIonisationModel();
      SetEmModel(ptr);
      ptr->SetLowEnergyLimit(0.0 * keV);
      ptr->SetHighEnergyLimit(3.0 * MeV);
    }

    AddEmModel(1, EmModel());

  } else if (name == "alpha") {

    if (!EmModel()) {
      auto ptr = new G4DNAQuadrupleIonisationModel();
      SetEmModel(ptr);
      ptr->SetLowEnergyLimit(0.0 * keV);
      ptr->SetHighEnergyLimit(23.0 * MeV);
    }

    AddEmModel(1, EmModel());

  } else if (name == "GenericIon") {

    // for carbon ions (12C6+)
    if (!EmModel()) {
      auto ptr = new G4DNAQuadrupleIonisationModel();
      SetEmModel(ptr);
      ptr->SetLowEnergyLimit(0.0 * keV);
      ptr->SetHighEnergyLimit(120.0 * MeV);
    }

    AddEmModel(1, EmModel());
  }

}

//------------------------------------------------------------------------------
void G4DNAQuadrupleIonisation::PrintInfo()
{
  if (EmModel(1)) {
    G4cout << " Total cross sections computed from " << EmModel(0)->GetName()
           << " and " << EmModel(1)->GetName() << " models" << G4endl;
  } else {
    G4cout << " Total cross sections computed from "
           << EmModel()->GetName() << G4endl;
  }
}
