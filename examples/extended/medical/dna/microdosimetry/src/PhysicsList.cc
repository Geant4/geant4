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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $ID$
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"

// Geant4-DNA MODELS

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAMillerGreenExcitationModel.hh"
#include "G4DNABornExcitationModel.hh"

#include "G4DNAIonisation.hh"
#include "G4DNABornIonisationModel.hh"
#include "G4DNARuddIonisationModel.hh"

#include "G4DNAChargeDecrease.hh"
#include "G4DNADingfelderChargeDecreaseModel.hh"

#include "G4DNAChargeIncrease.hh"
#include "G4DNADingfelderChargeIncreaseModel.hh"

#include "G4DNAAttachment.hh"
#include "G4DNAMeltonAttachmentModel.hh"

#include "G4DNAVibExcitation.hh"
#include "G4DNASancheExcitationModel.hh"

#include "G4DNAElectronSolvation.hh"

//

#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4VEmModel.hh"
#include "G4DummyModel.hh"
#include "G4eIonisation.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4BraggIonGasModel.hh"
#include "G4BetheBlochIonGasModel.hh"
#include "G4UrbanMscModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"

#include "G4ElectronCapture.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1*micrometer;
  fCutForGamma     = defaultCutValue;
  fCutForElectron  = defaultCutValue;
  fCutForPositron  = defaultCutValue;
  fCutForProton    = defaultCutValue;

  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructBosons()
{ 
  // gamma
  G4Gamma::GammaDefinition();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//DNA
#include "G4DNAGenericIonsManager.hh"
//ENDDNA

void PhysicsList::ConstructBarions()
{
  //  baryons
  G4Proton::ProtonDefinition();
  G4GenericIon::GenericIonDefinition();

  // Geant4 DNA new particles
  G4DNAGenericIonsManager * genericIonsManager;
  genericIonsManager=G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructEM()
{

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();

  while( (*particleIterator)() )
  {

    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    // *********************************
    // 1) Processes for the World region
    // *********************************

    if (particleName == "e-") {

      // STANDARD msc is active in the world
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->AddEmModel(1,new G4UrbanMscModel());
      pmanager->AddProcess(msc, -1, 1, -1);

      // STANDARD ionisation is active in the world
      G4eIonisation* eion = new G4eIonisation();
      eion->SetEmModel(new G4MollerBhabhaModel(), 1);
      pmanager->AddProcess(eion, -1, 2, 2);

      // DNA elastic is not active in the world 
      G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DummyModel(),1);
      pmanager->AddDiscreteProcess(theDNAElasticProcess);

      // DNA excitation is not active in the world 
      G4DNAExcitation* dnaex = new G4DNAExcitation("e-_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel(),1);
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world 
      G4DNAIonisation* dnaioni = new G4DNAIonisation("e-_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel(),1); 
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA attachment is not active in the world 
      G4DNAAttachment* dnaatt = new G4DNAAttachment("e-_G4DNAAttachment");
      dnaatt->SetEmModel(new G4DummyModel(),1); 
      pmanager->AddDiscreteProcess(dnaatt);

      // DNA vib. excitation is not active in the world 
      G4DNAVibExcitation* dnavib =
          new G4DNAVibExcitation("e-_G4DNAVibExcitation");
      dnavib->SetEmModel(new G4DummyModel(),1); 
      pmanager->AddDiscreteProcess(dnavib);

      // THE FOLLOWING PROCESS WILL KILL ALL ELECTRONS BELOW A
      // SELECTED ENERY THRESHOLD
      // Capture of low-energy e-
      G4ElectronCapture* ecap = new G4ElectronCapture("Target", 7.4*eV);
      pmanager->AddDiscreteProcess(ecap);
      // 7.4 eV is compatible with the validity range of the Champion's model
      
    } else if ( particleName == "proton" ) {

      // STANDARD msc is active in the world 
      G4hMultipleScattering* msc = new G4hMultipleScattering();
      msc->AddEmModel(1,new G4UrbanMscModel());
      pmanager->AddProcess(msc, -1, 1, -1);

      // STANDARD ionisation is active in the world 
      G4hIonisation* hion = new G4hIonisation();
      hion->SetEmModel(new G4BraggIonGasModel(), 1);
      hion->SetEmModel(new G4BetheBlochIonGasModel(), 2);
      pmanager->AddProcess(hion, -1, 2, 2);

      // DNA excitation is not active in the world 
      G4DNAExcitation* dnaex = new G4DNAExcitation("proton_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel(),1);
      dnaex->SetEmModel(new G4DummyModel(),2);
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world 
      G4DNAIonisation* dnaioni = new G4DNAIonisation("proton_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel(),1); 
      dnaioni->SetEmModel(new G4DummyModel(),2); 
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA charge decrease is ACTIVE in the world since
      // no corresponding STANDARD process exist
      pmanager->AddDiscreteProcess(
          new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"));

    } else if ( particleName == "hydrogen" ) {

      // DNA processes are ACTIVE in the world since
      // no corresponding STANDARD processes exist
      pmanager->AddDiscreteProcess(
          new G4DNAIonisation("hydrogen_G4DNAIonisation"));
      pmanager->AddDiscreteProcess(
          new G4DNAExcitation("hydrogen_G4DNAExcitation"));
      pmanager->AddDiscreteProcess(
          new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease"));

    } else if (particleName == "GenericIon") { 
    
      // THIS IS NEEDED FOR STANDARD ALPHA G4ionIonisation PROCESS

      // STANDARD msc is active in the world 
      G4hMultipleScattering* msc = new G4hMultipleScattering();
      msc->AddEmModel(1, new G4UrbanMscModel());
      pmanager->AddProcess(msc, -1, 1, -1);

      // STANDARD ionisation is active in the world 
      G4ionIonisation* hion = new G4ionIonisation();
      hion->SetEmModel(new G4BraggIonGasModel(),1);
      hion->SetEmModel(new G4BetheBlochIonGasModel(), 2);
      pmanager->AddProcess(hion, -1, 2, 2);

    } else if ( particleName == "alpha" ) {

      // STANDARD msc is active in the world 
      G4hMultipleScattering* msc = new G4hMultipleScattering();
      msc->AddEmModel(1, new G4UrbanMscModel());
      pmanager->AddProcess(msc, -1, 1, -1);

      // STANDARD ionisation is active in the world 
      G4ionIonisation* hion = new G4ionIonisation();
      hion->SetEmModel(new G4BraggIonGasModel(),1);
      hion->SetEmModel(new G4BetheBlochIonGasModel(), 2);
      pmanager->AddProcess(hion, -1, 2, 2);

      // DNA excitation is not active in the world 
      G4DNAExcitation* dnaex = new G4DNAExcitation("alpha_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel(),1);
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world 
      G4DNAIonisation* dnaioni = new G4DNAIonisation("alpha_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel(),1); 
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA charge decrease is ACTIVE in the world since no
      // corresponding STANDARD process exist
      pmanager->AddDiscreteProcess(
          new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease"));

    } else if ( particleName == "alpha+" ) {

      // DNA processes are ACTIVE in the world since no
      // corresponding STANDARD processes exist
      pmanager->AddDiscreteProcess(
          new G4DNAExcitation("alpha+_G4DNAExcitation"));
      pmanager->AddDiscreteProcess(
          new G4DNAIonisation("alpha+_G4DNAIonisation"));
      pmanager->AddDiscreteProcess(
          new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease"));
      pmanager->AddDiscreteProcess(
          new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease"));

    } else if ( particleName == "helium" ) {

      // DNA processes are ACTIVE in the world since no
      // corresponding STANDARD processes exist
      pmanager->AddDiscreteProcess(
          new G4DNAExcitation("helium_G4DNAExcitation"));
      pmanager->AddDiscreteProcess(
          new G4DNAIonisation("helium_G4DNAIonisation"));
      pmanager->AddDiscreteProcess(
          new G4DNAChargeIncrease("helium_G4DNAChargeIncrease"));

    }
  }

  // **************************************
  // 2) Define processes for Target region 
  // **************************************

  // STANDARD EM processes should be inactivated when
  // corresponding DNA processes are used
  // - STANDARD EM e- processes are inactivated below 1 MeV
  // - STANDARD EM proton & alpha processes are inactivated below
  //   standEnergyLimit
  G4double standEnergyLimit = 9.9*MeV;
  //

  G4double massFactor = 1.0079/4.0026;
  G4EmConfigurator* em_config =
      G4LossTableManager::Instance()->EmConfigurator();

  G4VEmModel* mod;

  // *** e-

  // ---> STANDARD EM processes are inactivated below 1 MeV

  mod =  new G4UrbanMscModel();
  mod->SetActivationLowEnergyLimit(1.*MeV);
  em_config->SetExtraEmModel("e-","msc",mod,"Target");

  mod = new G4MollerBhabhaModel();
  mod->SetActivationLowEnergyLimit(1*MeV);
  em_config->SetExtraEmModel("e-",
                             "eIoni",
                             mod,
                             "Target",
                             0.0,
                             100*TeV,
                             new G4UniversalFluctuation());

  // ---> DNA processes activated

  mod = new G4DNAChampionElasticModel();
  em_config->SetExtraEmModel("e-","e-_G4DNAElastic",
                             mod,"Target",7.4*eV,1.*MeV);

  mod = new G4DNABornIonisationModel();
  em_config->SetExtraEmModel("e-","e-_G4DNAIonisation",
                             mod,"Target",11.*eV,1.*MeV);
  // Note: valid from 11 eV to 0.999.. MeV then switch to std models at
  // higher energies ; same for other models

  mod = new G4DNABornExcitationModel();
  em_config->SetExtraEmModel("e-","e-_G4DNAExcitation",
                             mod,"Target",9.*eV,1.*MeV);

  mod = new G4DNAMeltonAttachmentModel();
  em_config->SetExtraEmModel("e-","e-_G4DNAAttachment",
                             mod,"Target",4.*eV,13.*eV);

  mod = new G4DNASancheExcitationModel();
  em_config->SetExtraEmModel("e-","e-_G4DNAVibExcitation",
                             mod,"Target",2.*eV,100.*eV);
  
  // *** proton

  // ---> STANDARD EM processes inactivated below standEnergyLimit

  // STANDARD msc is still active
  // Inactivate following STANDARD processes 

  mod = new G4BraggIonGasModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("proton","hIoni",
                             mod,"Target",0.0,2*MeV, new G4IonFluctuations());

  mod = new G4BetheBlochIonGasModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("proton","hIoni",
                             mod,"Target",2*MeV,100*TeV,
                             new G4UniversalFluctuation());

  // ---> DNA processes activated

  mod = new G4DNARuddIonisationModel(); 
  em_config->SetExtraEmModel("proton","proton_G4DNAIonisation",
                             mod,"Target",0.0,0.5*MeV);

  mod = new G4DNABornIonisationModel();
  em_config->SetExtraEmModel("proton","proton_G4DNAIonisation",
                             mod,"Target",0.5*MeV,10*MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("proton","proton_G4DNAExcitation",
                             mod,"Target",10*eV,0.5*MeV);

  mod = new G4DNABornExcitationModel();
  em_config->SetExtraEmModel("proton","proton_G4DNAExcitation",
                             mod,"Target",0.5*MeV,10*MeV);

  // *** alpha

  // ---> STANDARD EM processes inactivated below standEnergyLimit

  // STANDARD msc is still active
  // Inactivate following STANDARD processes 

  mod = new G4BraggIonGasModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha","ionIoni",
                             mod,"Target",0.0,2*MeV/massFactor,
                             new G4IonFluctuations());

  mod = new G4BetheBlochIonGasModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha","ionIoni",
                             mod,"Target",2*MeV/massFactor,100*TeV,
                             new G4UniversalFluctuation());

  // ---> DNA processes activated

  mod = new G4DNARuddIonisationModel(); 
  em_config->SetExtraEmModel("alpha","alpha_G4DNAIonisation",
                             mod,"Target",0.0,10*MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("alpha","alpha_G4DNAExcitation",
                             mod,"Target",1*keV,10*MeV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructGeneral()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  if (verboseLevel >0)
  {
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(fCutForGamma, "gamma");
  SetCutValue(fCutForElectron, "e-");
  SetCutValue(fCutForPositron, "e+");
  SetCutValue(fCutForProton, "proton");

  if (verboseLevel>0) { DumpCutValuesTable(); }
}
