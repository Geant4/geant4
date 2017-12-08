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
// $Id: G4MuIonisation.cc 107056 2017-11-01 14:52:32Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4MuIonisation
//
// Author:        Laszlo Urban
//
// Creation date: 30.09.1997
//
// Modifications:
//
// 08-04-98 remove 'tracking cut' of the ionizing particle (mma)
// 26-10-98 new stuff from R.Kokoulin + cleanup , L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 23-03-01 R.Kokoulin's correction is commented out, L.Urban
// 29-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 28-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 26-09-01 completion of RetrievePhysicsTable (mma)
// 29-10-01 all static functions no more inlined (mma)
// 07-11-01 correction(Tmax+xsection computation) L.Urban
// 08-11-01 particleMass becomes a local variable (mma)
// 10-05-02 V.Ivanchenko update to new design
// 04-12-02 V.Ivanchenko the low energy limit for Kokoulin model to 10 GeV
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 23-05-03 Define default integral + BohrFluctuations (V.Ivanchenko)
// 03-06-03 Add SetIntegral method to choose fluctuation model (V.Ivanchenko)
// 03-06-03 Fix initialisation problem for STD ionisation (V.Ivanchenko)
// 04-08-03 Set integral=false to be default (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 10-02-04 Calculation of radiative corrections using R.Kokoulin model (V.I.)
// 27-05-04 Set integral to be a default regime (V.Ivanchenko)
// 17-08-04 Utilise mu+ tables for mu- (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 12-08-05 SetStepLimits(0.2, 0.1*mm) (mma)
// 02-09-05 SetStepLimits(0.2, 1*mm) (V.Ivantchenko)
// 12-08-05 SetStepLimits(0.2, 0.1*mm) + integral off (V.Ivantchenko)
// 10-01-06 SetStepLimits -> SetStepFunction (V.Ivantchenko)
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MuIonisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4MuBetheBlochModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4IonFluctuations.hh"
#include "G4BohrFluctuations.hh"
#include "G4UnitsTable.hh"
#include "G4ICRU73QOModel.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4MuIonisation::G4MuIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(nullptr),
    theBaseParticle(nullptr),
    isInitialised(false)
{
  mass = ratio = 0;
  SetProcessSubType(fIonisation);
  SetSecondaryParticle(G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuIonisation::~G4MuIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && p.GetPDGMass() > 10.0*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuIonisation::MinPrimaryEnergy(const G4ParticleDefinition*,
					  const G4Material*,
					  G4double cut)
{
  G4double x = 0.5*cut/electron_mass_c2;
  G4double gam = x*ratio + std::sqrt((1. + x)*(1. + x*ratio*ratio));
  return mass*(gam - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuIonisation::InitialiseEnergyLossProcess(const G4ParticleDefinition* part,
                                                 const G4ParticleDefinition* bpart)
{
  if(!isInitialised) {

    theParticle = part;
    theBaseParticle = bpart;

    mass = theParticle->GetPDGMass();
    G4double q = theParticle->GetPDGCharge();

    G4EmParameters* param = G4EmParameters::Instance();
    G4double elow = 0.2*MeV;
    G4double emax = param->MaxKinEnergy();
    G4double ehigh = std::min(1*GeV, emax);

    // Bragg peak model
    if (!EmModel(0)) {
      if(q > 0.0) { SetEmModel(new G4BraggModel()); }
      else        { SetEmModel(new G4ICRU73QOModel()); }
    }
    EmModel(0)->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel(0)->SetHighEnergyLimit(elow); 
    AddEmModel(1, EmModel(0), new G4IonFluctuations());

    // high energy fluctuation model
    if (!FluctModel()) { SetFluctModel(new G4UniversalFluctuation()); }

    // moderate energy model
    if (!EmModel(1)) { SetEmModel(new G4BetheBlochModel()); }
    EmModel(1)->SetLowEnergyLimit(elow);
    EmModel(1)->SetHighEnergyLimit(ehigh);
    AddEmModel(2, EmModel(1), FluctModel());

    // high energy model
    if(ehigh < emax) {
      if (!EmModel(2)) { SetEmModel(new G4MuBetheBlochModel()); }
      EmModel(2)->SetLowEnergyLimit(ehigh);
      EmModel(2)->SetHighEnergyLimit(emax);
      AddEmModel(3, EmModel(2), FluctModel());
    }
    ratio = electron_mass_c2/mass;
    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuIonisation::PrintInfo()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuIonisation::ProcessDescription(std::ostream& out) const
{
  out << "<strong>Muon ionisation</strong>";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
