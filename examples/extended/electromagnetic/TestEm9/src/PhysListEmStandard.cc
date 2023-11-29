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
//---------------------------------------------------------------------------
//
// ClassName:   PhysListEmStandard
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "PhysListEmStandard.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4EmParameters.hh"
#include "G4EmBuilder.hh"
#include "G4LossTableManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4UrbanMscModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4GenericIon.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::PhysListEmStandard(G4int verbose)
  : G4VPhysicsConstructor("local")
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(verbose);
  SetPhysicsType(bElectromagnetic);
  SetVerboseLevel(verbose);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::~PhysListEmStandard()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructParticle()
{
  // minimal set of particles for EM physics
  G4EmBuilder::ConstructMinimalEmSet();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructProcess()
{
  if(verboseLevel > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4EmBuilder::PrepareEMPhysics();

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // processes used by several particles
  G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");
  G4NuclearStopping* pnuc(nullptr);

  // high energy limit for e+- scattering models and bremsstrahlung
  G4double highEnergyLimit = G4EmParameters::Instance()->MscEnergyLimit();

  // Add gamma EM Processes
  G4ParticleDefinition* particle = G4Gamma::Gamma();

  G4PhotoElectricEffect* pee = new G4PhotoElectricEffect();
  pee->SetEmModel(new G4LivermorePhotoElectricModel());
  ph->RegisterProcess(pee, particle);
  ph->RegisterProcess(new G4ComptonScattering(), particle);
  ph->RegisterProcess(new G4GammaConversion(), particle);
  ph->RegisterProcess(new G4RayleighScattering(), particle);

  // e-
  particle = G4Electron::Electron();

  G4eMultipleScattering* msc = new G4eMultipleScattering;
  G4UrbanMscModel* msc1 = new G4UrbanMscModel();
  G4WentzelVIModel* msc2 = new G4WentzelVIModel();
  msc1->SetHighEnergyLimit(highEnergyLimit);
  msc2->SetLowEnergyLimit(highEnergyLimit);
  msc->SetEmModel(msc1);
  msc->SetEmModel(msc2);

  G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel(); 
  G4CoulombScattering* ss = new G4CoulombScattering();
  ss->SetEmModel(ssm); 
  ss->SetMinKinEnergy(highEnergyLimit);
  ssm->SetLowEnergyLimit(highEnergyLimit);
  ssm->SetActivationLowEnergyLimit(highEnergyLimit);

  ph->RegisterProcess(msc, particle);
  ph->RegisterProcess(new G4eIonisation(), particle);
  ph->RegisterProcess(new G4eBremsstrahlung(), particle);
  ph->RegisterProcess(ss, particle);

  // e+
  particle = G4Positron::Positron();

  msc = new G4eMultipleScattering;
  msc1 = new G4UrbanMscModel();
  msc2 = new G4WentzelVIModel();
  msc1->SetHighEnergyLimit(highEnergyLimit);
  msc2->SetLowEnergyLimit(highEnergyLimit);
  msc->SetEmModel(msc1);
  msc->SetEmModel(msc2);

  ssm = new G4eCoulombScatteringModel(); 
  ss = new G4CoulombScattering();
  ss->SetEmModel(ssm); 
  ss->SetMinKinEnergy(highEnergyLimit);
  ssm->SetLowEnergyLimit(highEnergyLimit);
  ssm->SetActivationLowEnergyLimit(highEnergyLimit);

  ph->RegisterProcess(msc, particle);
  ph->RegisterProcess(new G4eIonisation(), particle);
  ph->RegisterProcess(new G4eBremsstrahlung(), particle);
  ph->RegisterProcess(new G4eplusAnnihilation(), particle);
  ph->RegisterProcess(ss, particle);

  // generic ion
  particle = G4GenericIon::GenericIon();
  G4ionIonisation* ionIoni = new G4ionIonisation();
  ph->RegisterProcess(hmsc, particle);
  ph->RegisterProcess(ionIoni, particle);

  // muons, hadrons ions
  G4EmBuilder::ConstructCharged(hmsc, pnuc);

  // extra configuration
  G4EmModelActivator mact(GetPhysicsName());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
