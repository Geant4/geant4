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
//------------------ G4PolarizedPhotoElectric physics process --

#include "G4PolarizedPhotoElectric.hh"

#include "G4PolarizedPhotoElectricModel.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4EmParameters.hh"

G4PolarizedPhotoElectric::G4PolarizedPhotoElectric(const G4String& processName,
                                                   G4ProcessType type)
  : G4VEmProcess(processName, type)
  , fIsInitialised(false)
{
  SetBuildTableFlag(false);
  SetSecondaryParticle(G4Electron::Electron());
  SetProcessSubType(fPhotoElectricEffect);
}

G4PolarizedPhotoElectric::~G4PolarizedPhotoElectric() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedPhotoElectric::ProcessDescription(std::ostream& out) const
{
  out << "Polarized model for photo-electric effect.\n";

  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4PolarizedPhotoElectric::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedPhotoElectric::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!fIsInitialised)
  {
    fIsInitialised = true;
    if(!EmModel())
    {
      SetEmModel(new G4PolarizedPhotoElectricModel);
    }
    G4EmParameters* param = G4EmParameters::Instance();
    EmModel()->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel()->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, EmModel());
  }
}
