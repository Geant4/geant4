//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4hhIonisation.cc,v 1.2 2005/11/29 08:13:48 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4hhIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 30.09.2005
//
// Modifications:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hhIonisation.hh"
#include "G4BraggNoDeltaModel.hh"
#include "G4BetheBlochNoDeltaModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4BohrFluctuations.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hhIonisation::G4hhIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(0),
    theBaseParticle(0),
    isInitialised(false)
{
  minKinEnergy = 0.1*keV;
  SetDEDXBinning(120);
  SetLambdaBinning(120);
  SetMinKinEnergy(minKinEnergy);
  SetMaxKinEnergy(100.0*TeV);
  SetVerboseLevel(2);
  mass = 0.0;
  ratio = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hhIonisation::~G4hhIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hhIonisation::InitialiseEnergyLossProcess(const G4ParticleDefinition* part,
                                                 const G4ParticleDefinition* bpart)
{
  if(isInitialised) return;

  theParticle = part;
  if(bpart) G4cout << "G4hhIonisation::InitialiseEnergyLossProcess WARNING: no "
                   << "base particle should be defined for the process "
		   << GetProcessName() << G4endl;

  SetBaseParticle(0);
  SetSecondaryParticle(G4Electron::Electron());
  mass  = theParticle->GetPDGMass();
  ratio = electron_mass_c2/mass;
  eth = 2.0*MeV*mass/proton_mass_c2;
  flucModel = new G4BohrFluctuations();

  G4int nm = 1;

  if(eth > minKinEnergy) {
    G4VEmModel* em = new G4BraggNoDeltaModel();
    em->SetLowEnergyLimit(minKinEnergy);
    em->SetHighEnergyLimit(eth);
    AddEmModel(nm, em, flucModel);
    nm++;
  }

  G4VEmModel* em1 = new G4BetheBlochNoDeltaModel();
  em1->SetLowEnergyLimit(std::max(eth,minKinEnergy));
  em1->SetHighEnergyLimit(100.0*TeV);
  AddEmModel(nm, em1, flucModel);

  SetStepLimits(0.1, 0.1*mm);

  G4cout << "G4hhIonisation is initialised: nm= " << nm << G4endl;

  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hhIonisation::PrintInfo()
{
  G4cout << "      Delta-ray will not be produced; "
         << "Bether-Bloch model for E > " << std::max(eth,minKinEnergy)
	 << G4endl;
  if(eth > minKinEnergy) G4cout
	 << "      ICRU49 parametrisation scaled from protons below.";
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
