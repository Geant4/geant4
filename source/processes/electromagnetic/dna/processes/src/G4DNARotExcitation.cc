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
#include "G4DNARotExcitation.hh"
#include "G4LEPTSRotExcitationModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4DNARotExcitation::G4DNARotExcitation(const G4String& processName) 
  : G4VEmProcess( processName ),isInitialised(false)
{
  SetBuildTableFlag(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4DNARotExcitation::IsApplicable(const G4ParticleDefinition& particleDef)
{
  return( &particleDef == G4Electron::Electron() ||
	  &particleDef == G4Positron::Positron() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4DNARotExcitation::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) 
  {
    isInitialised = true;
    if(nullptr == EmModel(0)) SetEmModel(new G4LEPTSRotExcitationModel);
    EmModel(0)->SetLowEnergyLimit(0.1*CLHEP::eV);
    EmModel(0)->SetHighEnergyLimit(15.*CLHEP::MeV);
    AddEmModel(1, EmModel(0));
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNARotExcitation::ProcessDescription(std::ostream& out) const
{
  out << "  DNA Rotation Excitation";
  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

