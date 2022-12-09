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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4MuonToMuonPairProduction
//
// Author:        Siddharth Yajaman on the base of Vladamir Ivantchenko/Laszlo Urban code
//
// Creation date: 12.07.2022
//
// Modifications:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MuonToMuonPairProduction.hh"
#include "G4SystemOfUnits.hh"
#include "G4MuonPlus.hh"
#include "G4VEmModel.hh"
#include "G4MuonToMuonPairProductionModel.hh"
#include "G4ElementData.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonToMuonPairProduction::G4MuonToMuonPairProduction(const G4String& name)
  : G4MuPairProduction(name)
{
  SetProcessSubType(49);
  SetSecondaryParticle(G4MuonPlus::MuonPlus());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonToMuonPairProduction::InitialiseEnergyLossProcess(
                         const G4ParticleDefinition* part,
			 const G4ParticleDefinition*)
{
  if (isInitialised) { return; }
  isInitialised = true;

  theParticle = part;
  lowestKinEnergy = std::max(lowestKinEnergy, part->GetPDGMass()*8.0);

  G4VEmModel* mod = EmModel(0);
  if(nullptr == mod) {
    mod = new G4MuonToMuonPairProductionModel(part); 
    SetEmModel(mod);
  }

  G4VEmFluctuationModel* fm = nullptr;
  G4EmParameters* param = G4EmParameters::Instance();
  mod->SetLowEnergyLimit(param->MinKinEnergy());
  mod->SetHighEnergyLimit(param->MaxKinEnergy());
  mod->SetSecondaryThreshold(param->MuHadBremsstrahlungTh());
  AddEmModel(1, mod, fm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonToMuonPairProduction::ProcessDescription(std::ostream& out) const
{
  out << "  Muon-anti-muon pair production by muons";
  G4MuPairProduction::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
