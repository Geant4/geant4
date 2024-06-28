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

#include "G4DNAChargeDecrease.hh"
#include "G4SystemOfUnits.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4DNADingfelderChargeDecreaseModel.hh"
#include "G4DNAIonChargeDecreaseModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAChargeDecrease::G4DNAChargeDecrease(const G4String& processName,
                                         G4ProcessType type) :
    G4VEmProcess(processName, type) 
{
  SetBuildTableFlag(false);
  SetProcessSubType(fLowEnergyChargeDecrease);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNAChargeDecrease::IsApplicable(const G4ParticleDefinition& p)
{
  G4String name = p.GetParticleName();
  return (p.IsGeneralIon() ||
	  name == "proton" || name == "alpha+" || name == "alpha");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAChargeDecrease::InitialiseProcess(const G4ParticleDefinition* p)
{
  if (!isInitialised) {
    isInitialised = true;

    auto mod = EmModel();
    if (nullptr == mod) { 
    
      G4String name = p->GetParticleName();

      if (name == "proton" || name == "alpha+" || name == "alpha") {
        SetEmModel(new G4DNADingfelderChargeDecreaseModel());
      } else {
        SetEmModel(new G4DNAIonChargeDecreaseModel());
      }
    }
    G4EmParameters* param = G4EmParameters::Instance();
    EmModel(0)->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel(0)->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, EmModel(0));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAChargeDecrease::ProcessDescription(std::ostream& out) const
{
  out << " Total cross sections computed from " << EmModel()->GetName()
      << " model" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
