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
// $Id: G4alphaIonisation.cc 107058 2017-11-01 14:54:12Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4alphaIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 28.10.2009 created from G4ionIonisation
//
// Modifications:
//
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4alphaIonisation.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Alpha.hh"
#include "G4BraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4alphaIonisation::G4alphaIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(nullptr),
    isInitialised(false)
{
  G4Exception("G4alphaIonisation::G4alphaIonisation","em0007",JustWarning,
	      " The process is not ready for use - incorrect results are expected");
  SetLinearLossLimit(0.02);
  SetProcessSubType(fIonisation);
  mass = 0.0;
  ratio = 0.0;
  eth = 8*MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4alphaIonisation::~G4alphaIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4alphaIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (!p.IsShortLived() &&
	  std::abs(p.GetPDGCharge()/CLHEP::eplus - 2) < 0.01);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4alphaIonisation::MinPrimaryEnergy(const G4ParticleDefinition*, 
					     const G4Material*, 
					     G4double cut)
{
  G4double x = 0.5*cut/electron_mass_c2;
  G4double gam = x*ratio + std::sqrt((1. + x)*(1. + x*ratio*ratio));
  return mass*(gam - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4alphaIonisation::InitialiseEnergyLossProcess(
		      const G4ParticleDefinition* part,
		      const G4ParticleDefinition* bpart)
{
  if(!isInitialised) {

    theParticle = part;
    G4String pname = part->GetParticleName();

    // define base particle
    const G4ParticleDefinition* theBaseParticle = nullptr;
    if(bpart == 0) { 
      if(pname != "alpha") { theBaseParticle = G4Alpha::Alpha(); }
    } else { theBaseParticle = bpart; }

    mass  = part->GetPDGMass();
    ratio = mass/electron_mass_c2/mass;

    SetBaseParticle(theBaseParticle);
    SetSecondaryParticle(G4Electron::Electron());

    if (!EmModel(0)) { SetEmModel(new G4BraggIonModel()); }

    G4EmParameters* param = G4EmParameters::Instance();
    G4double emin = param->MinKinEnergy();
    EmModel(0)->SetLowEnergyLimit(emin);

    // model limit defined for alpha
    eth = (EmModel(0)->HighEnergyLimit())*ratio;
    EmModel(0)->SetHighEnergyLimit(eth);
    AddEmModel(1, EmModel(0), new G4IonFluctuations());

    if (!FluctModel()) { SetFluctModel(new G4UniversalFluctuation()); }

    if (!EmModel(1)) { SetEmModel(new G4BetheBlochModel()); }  
    EmModel(1)->SetLowEnergyLimit(eth);
    EmModel(1)->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(2, EmModel(1), FluctModel());    

    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4alphaIonisation::PrintInfo()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4alphaIonisation::ProcessDescription(std::ostream& out) const
{
  out << "<strong>Alpha ionisation</strong>";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
