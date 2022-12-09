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
// GEANT4 Class file
//
//
// File name:     G4hIonisation
//
// Author:        Laszlo Urban
//
// Creation date: 30.05.1997
//
// Modified by Laszlo Urban, Michel Maire and Vladimir Ivanchenko
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hIonisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4EmStandUtil.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4ICRU73QOModel.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hIonisation::G4hIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    isInitialised(false)
{
  SetProcessSubType(fIonisation);
  SetSecondaryParticle(G4Electron::Electron());
  mass = 0.0;
  ratio = 0.0;
  eth = 2*CLHEP::MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hIonisation::~G4hIonisation() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4hIonisation::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hIonisation::MinPrimaryEnergy(const G4ParticleDefinition*,
					 const G4Material*,
					 G4double cut)
{
  G4double x = 0.5*cut/electron_mass_c2;
  G4double gam = x*ratio + std::sqrt((1. + x)*(1. + x*ratio*ratio));
  return mass*(gam - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....  

void G4hIonisation::InitialiseEnergyLossProcess(
		    const G4ParticleDefinition* part,
		    const G4ParticleDefinition* bpart)
{
  if(!isInitialised) {

    const G4ParticleDefinition* theBaseParticle = nullptr;
    G4String pname = part->GetParticleName();
    G4double q = part->GetPDGCharge();

    //G4cout << " G4hIonisation::InitialiseEnergyLossProcess " << pname 
    //   << "  " << bpart << G4endl;

    // define base particle
    if(part == bpart) { 
      theBaseParticle = nullptr;
    } else if(nullptr != bpart) { 
      theBaseParticle = bpart;

    } else if(pname == "proton" || pname == "anti_proton" || 
	      pname == "pi+" || pname == "pi-" || 
	      pname == "kaon+" || pname == "kaon-" || 
	      pname == "GenericIon" || pname == "alpha") { 
      // no base particles
      theBaseParticle = nullptr;

    } else {
      // select base particle 
      if(part->GetPDGSpin() == 0.0) {
	if(q > 0.0) { theBaseParticle = G4KaonPlus::KaonPlus(); }
	else { theBaseParticle = G4KaonMinus::KaonMinus(); }
      } else {
	if(q > 0.0) { theBaseParticle = G4Proton::Proton(); } 
	else { theBaseParticle = G4AntiProton::AntiProton(); }
      }
    }
    SetBaseParticle(theBaseParticle);

    // model limit defined for protons
    mass  = part->GetPDGMass();
    ratio = electron_mass_c2/mass;
    eth   = 2.0*MeV*mass/proton_mass_c2;

    G4EmParameters* param = G4EmParameters::Instance();
    G4double emin = param->MinKinEnergy();
    G4double emax = param->MaxKinEnergy();

    // define model of energy loss fluctuations
    if (nullptr == FluctModel()) {
      G4bool ion = (pname == "GenericIon" || pname == "alpha"); 
      SetFluctModel(G4EmStandUtil::ModelOfFluctuations(ion));
    }

    if (nullptr == EmModel(0)) { 
      if(q > 0.0) { SetEmModel(new G4BraggModel()); }
      else        { SetEmModel(new G4ICRU73QOModel()); }
    }
    // to compute ranges correctly we have to use low-energy
    // model even if activation limit is high
    EmModel(0)->SetLowEnergyLimit(emin);

    // high energy limit may be eth or DBL_MAX
    G4double emax1 = (EmModel(0)->HighEnergyLimit() < emax) ? eth : emax;
    EmModel(0)->SetHighEnergyLimit(emax1);
    AddEmModel(1, EmModel(0), FluctModel());
    
    // second model is used if the first does not cover energy range
    if(emax1 < emax) {
      if (nullptr == EmModel(1)) { SetEmModel(new G4BetheBlochModel()); }
      EmModel(1)->SetLowEnergyLimit(emax1);

      // for extremely heavy particles upper limit of the model
      // should be increased
      emax = std::max(emax, eth*10); 
      EmModel(1)->SetHighEnergyLimit(emax);
      AddEmModel(2, EmModel(1), FluctModel());  
    }
    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hIonisation::ProcessDescription(std::ostream& out) const
{
  out << "  Ionisation";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
