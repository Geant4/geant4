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
// $Id: G4hIonisation.cc 102525 2017-02-08 11:23:26Z gcosmo $
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
// Modifications:
//
// corrected by L.Urban on 24/09/97
// several bugs corrected by L.Urban on 13/01/98
// 07-04-98 remove 'tracking cut' of the ionizing particle, mma
// 22-10-98 cleanup L.Urban
// 02-02-99 bugs fixed , L.Urban
// 29-07-99 correction in BuildLossTable for low energy, L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 10-08-00 V.Ivanchenko change BuildLambdaTable, in order to
//          simulate energy losses of ions; correction to
//          cross section for particles with spin 1 is inserted as well
// 28-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 14-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 29-08-01 PostStepDoIt: correction for spin 1/2 (instead of 1) (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 25-09-01 completion of RetrievePhysicsTable() (mma)
// 29-10-01 all static functions no more inlined
// 08-11-01 Charge renamed zparticle; added to the dedx
// 27-03-02 Bug fix in scaling of lambda table (V.Ivanchenko)
// 09-04-02 Update calculation of tables for GenericIons (V.Ivanchenko)
// 30-04-02 V.Ivanchenko update to new design
// 04-12-02 Add verbose level definition (VI)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 23-05-03 Define default integral + BohrFluctuations (V.Ivanchenko)
// 03-06-03 Fix initialisation problem for STD ionisation (V.Ivanchenko)
// 04-08-03 Set integral=false to be default (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 27-05-04 Set integral to be a default regime (V.Ivanchenko) 
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 24-03-05 Optimize internal interfaces (V.Ivantchenko)
// 12-08-05 SetStepLimits(0.2, 0.1*mm) (mma)
// 10-01-06 SetStepLimits -> SetStepFunction (V.Ivanchenko)
// 26-05-06 scale negative particles from pi- and pbar,
//          positive from pi+ and p (VI)
// 14-01-07 use SetEmModel() and SetFluctModel() from G4VEnergyLossProcess (mma)
// 12-09-08 Removed CorrectionsAlongStep (VI)
// 27-05-10 Added G4ICRU73QOModel for anti-protons (VI)
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
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"
#include "G4BohrFluctuations.hh"
#include "G4UnitsTable.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4ICRU73QOModel.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4hIonisation::G4hIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    isInitialised(false)
{
  SetProcessSubType(fIonisation);
  SetSecondaryParticle(G4Electron::Electron());
  mass = 0.0;
  ratio = 0.0;
  eth = 2*MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hIonisation::~G4hIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4hIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && p.GetPDGMass() > 10.0*MeV &&
	 !p.IsShortLived());
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

    // standard base particles
    if(part == bpart || pname == "proton" ||
       pname == "anti_proton" || 
       pname == "pi+" || pname == "pi-" || 
       pname == "kaon+" || pname == "kaon-" || pname == "GenericIon"
       || pname == "He3" || pname == "alpha") 
      { 
	theBaseParticle = nullptr;
      }
    // select base particle 
    else if(bpart == nullptr) {

      if(part->GetPDGSpin() == 0.0) {
	if(q > 0.0) { theBaseParticle = G4KaonPlus::KaonPlus(); }
	else { theBaseParticle = G4KaonMinus::KaonMinus(); }
      } else {
	if(q > 0.0) { theBaseParticle = G4Proton::Proton(); } 
	else { theBaseParticle = G4AntiProton::AntiProton(); }
      }

      // base particle defined by interface
    } else { 
      theBaseParticle = bpart;
    }
    SetBaseParticle(theBaseParticle);

    mass  = part->GetPDGMass();
    ratio = electron_mass_c2/mass;
    eth   = 2.0*MeV*mass/proton_mass_c2;

    G4EmParameters* param = G4EmParameters::Instance();
    G4double emin = std::min(param->MinKinEnergy(), 0.1*eth);
    G4double emax = std::max(param->MaxKinEnergy(), 100*eth);

    if(emin != param->MinKinEnergy() || emax != param->MaxKinEnergy()) {
      SetMinKinEnergy(emin);
      SetMaxKinEnergy(emax);
      G4int bin = G4lrint(param->NumberOfBinsPerDecade()*std::log10(emax/emin));
      SetDEDXBinning(bin);
    }

    if (!EmModel(1)) { 
      if(q > 0.0) { SetEmModel(new G4BraggModel(),1); }
      else { SetEmModel(new G4ICRU73QOModel(),1); }
    }
    EmModel(1)->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel(1)->SetHighEnergyLimit(eth);
    AddEmModel(1, EmModel(1), new G4IonFluctuations());

    if (!FluctModel()) { SetFluctModel(new G4UniversalFluctuation()); }

    if (!EmModel(2)) { SetEmModel(new G4BetheBlochModel(),2); }
    EmModel(2)->SetLowEnergyLimit(eth);
    EmModel(2)->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(2, EmModel(2), FluctModel());  

    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hIonisation::PrintInfo()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
