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
#include "G4DNADissociation.hh"
#include "G4LEPTSDissociationModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4DNADissociation::G4DNADissociation(const G4String& processName) :
    G4VEmProcess(processName) 
{
} // constructor

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4DNADissociation::~G4DNADissociation()
= default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4DNADissociation::IsApplicable(const G4ParticleDefinition& particleDef)
{
  return (&particleDef == G4Electron::Electron()
      || &particleDef == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4DNADissociation::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised)
  {
    isInitialised = true;
    SetBuildTableFlag(false);
    if(EmModel() == nullptr)
    {
      SetEmModel(new G4LEPTSDissociationModel);
      EmModel()->SetLowEnergyLimit(0.1 * CLHEP::eV);
      EmModel()->SetHighEnergyLimit(15. * CLHEP::MeV);
    }
    AddEmModel(1, EmModel());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4DNADissociation::PrintInfo()
{
  G4cout << " Total cross sections computed from " << EmModel()->GetName()
         << " model" << G4endl;
}

