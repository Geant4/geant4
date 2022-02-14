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
// -------------------------------------------------------------------
//
// Geant4 Class file
//
// File name:     G4PolarizedBremsstrahlung
//
// Author:        Karim Laihem

#include "G4PolarizedBremsstrahlung.hh"

#include "G4EmParameters.hh"
#include "G4Gamma.hh"
#include "G4PolarizedBremsstrahlungModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedBremsstrahlung::G4PolarizedBremsstrahlung(const G4String& name)
  : G4eBremsstrahlung(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedBremsstrahlung::~G4PolarizedBremsstrahlung() {}

void G4PolarizedBremsstrahlung::ProcessDescription(std::ostream& out) const
{
  out << "Polarized model for bremsstrahlung.\n";
  G4eBremsstrahlung::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedBremsstrahlung::InitialiseEnergyLossProcess(
  const G4ParticleDefinition*, const G4ParticleDefinition*)
{
  if(!isInitialised)
  {
    isInitialised = true;
    G4VEmFluctuationModel* fm = nullptr;
    G4VEmModel* em        = new G4PolarizedBremsstrahlungModel;
    G4EmParameters* param = G4EmParameters::Instance();
    em->SetLowEnergyLimit(param->MinKinEnergy());
    em->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, em, fm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
