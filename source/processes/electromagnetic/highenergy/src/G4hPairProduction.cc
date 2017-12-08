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
// $Id: G4hPairProduction.cc 107367 2017-11-09 10:56:09Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4hPairProduction
//
// Author:        Vladimir Ivanchenko on base of model for muons
//
// Creation date: 01.03.2008
//
// Modifications:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hPairProduction.hh"
#include "G4SystemOfUnits.hh"
#include "G4hPairProductionModel.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4hPairProduction::G4hPairProduction(const G4String& name)
  : G4MuPairProduction(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hPairProduction::~G4hPairProduction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4hPairProduction::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && p.GetPDGMass() > 110.0*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hPairProduction::InitialiseEnergyLossProcess(
                         const G4ParticleDefinition* part,
			 const G4ParticleDefinition*)
{
  if (!isInitialised) {
    isInitialised = true;

    theParticle = part;
    if (!EmModel()) { SetEmModel(new G4hPairProductionModel(part)); }

    G4double limit = part->GetPDGMass()*8.;
    if(limit > lowestKinEnergy) { lowestKinEnergy = limit; }

    G4VEmFluctuationModel* fm = nullptr;
    G4EmParameters* param = G4EmParameters::Instance();
    EmModel()->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel()->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, EmModel(), fm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hPairProduction::ProcessDescription(std::ostream& out) const
{
  out << "<strong>Hadron pair production</strong>";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
