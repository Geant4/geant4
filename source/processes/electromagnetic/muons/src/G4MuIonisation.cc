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
// $Id: G4MuIonisation.cc,v 1.54 2007/05/22 17:35:58 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
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
// 10-02-04 Calculation of radiative corrections using R.Kokoulin model (V.Ivanchenko)
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
#include "G4Electron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4MuBetheBlochModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4BohrFluctuations.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4MuIonisation::G4MuIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(0),
    theBaseParticle(0),
    isInitialised(false)
{
  SetStepFunction(0.2, 1*mm);
  SetIntegral(true);
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuIonisation::~G4MuIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuIonisation::InitialiseEnergyLossProcess(const G4ParticleDefinition* part,
                                                 const G4ParticleDefinition* bpart)
{
  if(!isInitialised) {

    theParticle = part;
    theBaseParticle = bpart;

    mass = theParticle->GetPDGMass();
    SetSecondaryParticle(G4Electron::Electron());

    flucModel = new G4UniversalFluctuation();

    G4VEmModel* em = new G4BraggModel();
    em->SetLowEnergyLimit(0.1*keV);
    em->SetHighEnergyLimit(0.2*MeV);
    AddEmModel(1, em, flucModel);
    G4VEmModel* em1 = new G4BetheBlochModel();
    em1->SetLowEnergyLimit(0.2*MeV);
    em1->SetHighEnergyLimit(1.0*GeV);
    AddEmModel(2, em1, flucModel);
    G4VEmModel* em2 = new G4MuBetheBlochModel();
    em2->SetLowEnergyLimit(1.0*GeV);
    em2->SetHighEnergyLimit(100.0*TeV);
    AddEmModel(3, em2, flucModel);

    ratio = electron_mass_c2/mass;
    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuIonisation::PrintInfo()
{
  G4cout << "      Bether-Bloch model for E > 0.2 MeV, "
         << "parametrisation of Bragg peak below, "
         << G4endl;
  G4cout << "      radiative corrections for E > 1 GeV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




