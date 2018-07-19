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
// $Id: G4eIonisation.cc 107058 2017-11-01 14:54:12Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eIonisation
//
// Author:        Laszlo Urban
//
// Creation date: 20.03.1997
//
// Modifications:
//
// 07-04-98 remove 'tracking cut' of the ionizing particle, mma
// 04-09-98 new methods SetBining() PrintInfo()
// 07-09-98 Cleanup
// 02-02-99 correction inDoIt , L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 28-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 09-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 13-08-01 new function ComputeRestrictedMeandEdx()  (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 21-09-01 completion of RetrievePhysicsTable() (mma)
// 29-10-01 all static functions no more inlined (mma)
// 07-11-01 particleMass and Charge become local variables
// 26-03-02 change access to cuts in BuildLossTables (V.Ivanchenko)
// 30-04-02 V.Ivanchenko update to new design
// 23-12-02 Change interface in order to move to cut per region (VI)
// 26-12-02 Secondary production moved to derived classes (VI)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 23-05-03 Define default integral + BohrFluctuations (V.Ivanchenko)
// 03-06-03 Fix initialisation problem for STD ionisation (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 12-08-05 SetStepLimits(0.2, 0.1*mm) (mma)
// 02-09-05 Return SetStepLimits(1, 1*mm) (V.Ivantchenko)
// 10-01-06 SetStepLimits -> SetStepFunction (V.Ivantchenko)
// 14-01-07 use SetEmModel() and SetFluctModel() from G4VEnergyLossProcess (mma)
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eIonisation.hh"
#include "G4Electron.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4BohrFluctuations.hh"
#include "G4UnitsTable.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eIonisation::G4eIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theElectron(G4Electron::Electron()),
    isElectron(true),
    isInitialised(false)
{
  SetProcessSubType(fIonisation);
  SetSecondaryParticle(theElectron);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eIonisation::~G4eIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eIonisation::MinPrimaryEnergy(const G4ParticleDefinition*,
					 const G4Material*,
					 G4double cut)
{
  G4double x = cut;
  if(isElectron) x += cut;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4eIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == theElectron || &p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisation::InitialiseEnergyLossProcess(
		    const G4ParticleDefinition* part,
		    const G4ParticleDefinition*)
{
  if(!isInitialised) {
    if(part != theElectron) { isElectron = false; }
    if (!EmModel(0)) { SetEmModel(new G4MollerBhabhaModel()); }
    G4EmParameters* param = G4EmParameters::Instance();
    EmModel(0)->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel(0)->SetHighEnergyLimit(param->MaxKinEnergy());
    if (!FluctModel()) { SetFluctModel(new G4UniversalFluctuation()); }
                
    AddEmModel(1, EmModel(), FluctModel());
    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisation::PrintInfo()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eIonisation::ProcessDescription(std::ostream& out) const
{
  out << "<strong>Ionisation</strong>";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
