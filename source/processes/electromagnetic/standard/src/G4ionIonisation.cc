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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ionIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.05.2002
//
// Modifications:
//
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 18-04-03 Use IonFluctuations (V.Ivanchenko)
// 03-08-03 Add effective charge (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 27-05-04 Set integral to be a default regime (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 10-01-06 SetStepLimits -> SetStepFunction (V.Ivantchenko)
// 10-05-06 Add a possibility to download user data (V.Ivantchenko)
// 13-05-06 Add data for light ion stopping in water (V.Ivantchenko)
// 14-01-07 use SetEmModel() and SetFluctModel() from G4VEnergyLossProcess (mma)
// 16-05-07 Add data for light ion stopping only for GenericIon (V.Ivantchenko)
// 07-11-07 Fill non-ionizing energy loss (V.Ivantchenko)
// 12-09-08 Removed InitialiseMassCharge and CorrectionsAlongStep (VI)
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ionIonisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4GenericIon.hh"
#include "G4BraggModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4LossTableManager.hh"
#include "G4WaterStopping.hh"
#include "G4EmCorrections.hh"
#include "G4EmParameters.hh"
#include "G4EmStandUtil.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionIonisation::G4ionIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(nullptr),
    isInitialised(false),
    stopDataActive(false)
{
  SetLinearLossLimit(0.02);
  SetProcessSubType(fIonisation);
  SetSecondaryParticle(G4Electron::Electron());
  corr = G4LossTableManager::Instance()->EmCorrections();
  eth = 2*CLHEP::MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionIonisation::~G4ionIonisation() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4ionIonisation::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionIonisation::MinPrimaryEnergy(const G4ParticleDefinition* p, 
					   const G4Material*, 
					   G4double cut)
{
  return p->GetPDGMass()*(std::sqrt(1. + 0.5*cut/CLHEP::electron_mass_c2) - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::InitialiseEnergyLossProcess(
		      const G4ParticleDefinition* part,
		      const G4ParticleDefinition* bpart)
{
  const G4ParticleDefinition* ion = G4GenericIon::GenericIon();

  if(!isInitialised) {
    theParticle = part;

    // define base particle
    const G4ParticleDefinition* theBaseParticle = nullptr;
    const G4int pdg = part->GetPDGEncoding();

    if(part == bpart) { 
      theBaseParticle = nullptr;
    } else if(nullptr != bpart) { 
      theBaseParticle = bpart;
    } else if(part == ion || pdg == 1000020040) { 
      theBaseParticle = nullptr;
    } else { 
      theBaseParticle = ion; 
    }
    SetBaseParticle(theBaseParticle);

    // model limit defined for protons
    eth = 2*CLHEP::MeV*part->GetPDGMass()/CLHEP::proton_mass_c2;

    G4EmParameters* param = G4EmParameters::Instance();
    G4double emin = param->MinKinEnergy();
    G4double emax = param->MaxKinEnergy();

    // define model of energy loss fluctuations
    if (nullptr == FluctModel()) {
      SetFluctModel(G4EmStandUtil::ModelOfFluctuations(true));
    }

    if (nullptr == EmModel(0)) { SetEmModel(new G4BraggIonModel()); }
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

      // Add ion stoping tables for Generic Ion if the default 
      // model is used (with eth ~= 2 MeV)
      if(part == ion && (EmModel(1)->GetName() == "BetheBloch" ||
			 EmModel(1)->GetName() == "BetheBlochGasIon")) {
	stopDataActive = true;
	G4WaterStopping  ws(corr, true);
	corr->SetIonisationModels(EmModel(0), EmModel(1));
      }
    }
    isInitialised = true;
  }
  // reinitialisation of corrections for the new run
  if(part == ion) { corr->InitialiseForNewRun(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::StreamProcessInfo(std::ostream& out) const
{
  if (stopDataActive && G4GenericIon::GenericIon() == theParticle) {
    out << "      Stopping Power data for " 
	<< corr->GetNumberOfStoppingVectors()
	<< " ion/material pairs" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::AddStoppingData(G4int Z, G4int A,
				      const G4String& mname,
				      G4PhysicsVector* dVector)
{
  corr->AddStoppingData(Z, A, mname, dVector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::ProcessDescription(std::ostream& out) const
{
  out << "  Ion ionisation";
  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 
