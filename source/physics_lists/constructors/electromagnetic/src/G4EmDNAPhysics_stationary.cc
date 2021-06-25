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

#include "G4EmDNAPhysics_stationary.hh"

#include "G4SystemOfUnits.hh"

#include "G4DNAGenericIonsManager.hh"

// *** Processes and models for Geant4-DNA

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAIonElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"

// particles

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4EmParameters.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_stationary);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4EmDNAPhysics_stationary::G4EmDNAPhysics_stationary(G4int ver, const G4String& nam)
  : G4EmDNAPhysics(ver, nam)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDNAStationary(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics_stationary::~G4EmDNAPhysics_stationary()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics_stationary::ConstructProcess()
{
  if(verboseLevel > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  ConstructGammaPositronProcesses();

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4DNAGenericIonsManager* genericIonsManager
    = G4DNAGenericIonsManager::Instance();

  // e-
  G4ParticleDefinition* part = G4Electron::Electron();

  // *** Elastic scattering ***      
  G4DNAElastic* theDNAElastic = new G4DNAElastic("e-_G4DNAElastic");
  theDNAElastic->SetEmModel(new G4DNAChampionElasticModel());
  ph->RegisterProcess(theDNAElastic, part);

  // *** Excitation ***
  G4DNAExcitation* theDNAExc = new G4DNAExcitation("e-_G4DNAExcitation");
  G4DNABornExcitationModel* modB = new G4DNABornExcitationModel();
  theDNAExc->SetEmModel(modB);
  modB->SelectStationary(true);
  ph->RegisterProcess(theDNAExc, part);

  // *** Ionisation ***
  G4DNAIonisation* theDNAIoni = new G4DNAIonisation("e-_G4DNAIonisation");
  G4DNABornIonisationModel* modI = new G4DNABornIonisationModel();
  theDNAIoni->SetEmModel(modI);
  modI->SelectStationary(true);
  ph->RegisterProcess(theDNAIoni, part);

  // *** Vibrational excitation ***
  G4DNAVibExcitation* theDNAVibExc = 
    new G4DNAVibExcitation("e-_G4DNAVibExcitation");
  G4DNASancheExcitationModel* modS = new G4DNASancheExcitationModel();
  theDNAVibExc->SetEmModel(modS);
  modS->SelectStationary(true);
  ph->RegisterProcess(theDNAVibExc, part);
          
  // *** Attachment ***
  G4DNAAttachment* theDNAAttach = new G4DNAAttachment("e-_G4DNAAttachment");
  G4DNAMeltonAttachmentModel* modM = new G4DNAMeltonAttachmentModel();
  theDNAAttach->SetEmModel(modM);
  modM->SelectStationary(true);
  ph->RegisterProcess(theDNAAttach, part);
          
  // proton
  part = G4Proton::Proton();
      
  // *** Elastic ***
  theDNAElastic = new G4DNAElastic("proton_G4DNAElastic");
  G4DNAIonElasticModel* modE = new G4DNAIonElasticModel();
  theDNAElastic->SetEmModel(modE);
  modE->SelectStationary(true);
  ph->RegisterProcess(theDNAElastic, part);
            
  // *** Excitation ***
  theDNAExc = new G4DNAExcitation("proton_G4DNAExcitation");
  G4DNAMillerGreenExcitationModel* modMGE = 
    new G4DNAMillerGreenExcitationModel();
  modMGE->SetLowEnergyLimit(10*eV);
  modMGE->SetHighEnergyLimit(500*keV);
  modMGE->SelectStationary(true);
  theDNAExc->SetEmModel(modMGE);
  G4DNABornExcitationModel* modBE = new G4DNABornExcitationModel();
  modBE->SetLowEnergyLimit(500*keV);
  modBE->SetHighEnergyLimit(100*MeV);
  modBE->SelectStationary(true);
  theDNAExc->SetEmModel(modBE);      
  ph->RegisterProcess(theDNAExc, part);

  // *** Ionisation ***
  theDNAIoni = new G4DNAIonisation("proton_G4DNAIonisation");
  G4DNARuddIonisationModel* modRI = new G4DNARuddIonisationModel();
  modRI->SetLowEnergyLimit(0*eV);
  modRI->SetHighEnergyLimit(500*keV);
  modRI->SelectStationary(true);
  theDNAIoni->SetEmModel(modRI);
  G4DNABornIonisationModel* modBI = new G4DNABornIonisationModel();
  modBI->SetLowEnergyLimit(500*keV);
  modBI->SetHighEnergyLimit(100*MeV);
  modBI->SelectStationary(true);
  theDNAIoni->SetEmModel(modBI);
  ph->RegisterProcess(theDNAIoni, part);

  // *** Charge decrease ***
  G4DNAChargeDecrease* theDNAChargeDecreaseProcess = 
    new G4DNAChargeDecrease("proton_G4DNAChargeDecrease");
  G4DNADingfelderChargeDecreaseModel* modDCD =
    new G4DNADingfelderChargeDecreaseModel();
  modDCD->SelectStationary(true);
  theDNAChargeDecreaseProcess->SetEmModel(modDCD);
  ph->RegisterProcess(theDNAChargeDecreaseProcess, part);

  // hydrogen
  part = genericIonsManager->GetIon("hydrogen");
      
  // *** Elastic ***
  theDNAElastic = new G4DNAElastic("hydrogen_G4DNAElastic");
  G4DNAIonElasticModel* modEI = new G4DNAIonElasticModel();
  modEI->SelectStationary(true);
  theDNAElastic->SetEmModel(modEI);
  ph->RegisterProcess(theDNAElastic, part);

  // *** Excitation ***
  theDNAExc = new G4DNAExcitation("hydrogen_G4DNAExcitation");
  modMGE = new G4DNAMillerGreenExcitationModel();
  modMGE->SelectStationary(true);
  theDNAExc->SetEmModel(modMGE);
  ph->RegisterProcess(theDNAExc, part);

  // *** Ionisation ***
  theDNAIoni = new G4DNAIonisation("hydrogen_G4DNAIonisation");
  modRI = new G4DNARuddIonisationModel();
  theDNAIoni->SetEmModel(new G4DNARuddIonisationModel());
  modRI->SelectStationary(true);
  ph->RegisterProcess(theDNAIoni, part);

  // *** Charge increase ***
  G4DNAChargeIncrease* theDNAChargeIncreaseProcess = 
    new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease");
  G4DNADingfelderChargeIncreaseModel* modDCI = 
    new G4DNADingfelderChargeIncreaseModel();
  modDCI->SelectStationary(true);
  theDNAChargeIncreaseProcess->SetEmModel(modDCI);
  ph->RegisterProcess(theDNAChargeIncreaseProcess, part);

  // alpha++
  part = G4Alpha::Alpha();
      
  // *** Elastic ***    
  theDNAElastic = new G4DNAElastic("alpha_G4DNAElastic");
  modEI = new G4DNAIonElasticModel();
  modEI->SelectStationary(true);
  theDNAElastic->SetEmModel(modEI);
  ph->RegisterProcess(theDNAElastic, part);
      
  // *** Excitation ***
  theDNAExc = new G4DNAExcitation("alpha_G4DNAExcitation");
  modMGE = new G4DNAMillerGreenExcitationModel();
  modMGE->SelectStationary(true);
  theDNAExc->SetEmModel(modMGE);
  ph->RegisterProcess(theDNAExc, part);
            
  // *** Ionisation ***
  theDNAIoni = new G4DNAIonisation("alpha_G4DNAIonisation");
  modRI = new G4DNARuddIonisationModel();
  theDNAIoni->SetEmModel(new G4DNARuddIonisationModel());
  modRI->SelectStationary(true);
  ph->RegisterProcess(theDNAIoni, part);

  // *** Charge decrease ***
  theDNAChargeDecreaseProcess = 
    new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease");
  modDCD = new G4DNADingfelderChargeDecreaseModel();
  modDCD->SelectStationary(true);
  theDNAChargeDecreaseProcess->SetEmModel(modDCD);
  ph->RegisterProcess(theDNAChargeDecreaseProcess, part);

  // alpha+
  part = genericIonsManager->GetIon("alpha+");
      
  // *** Elastic ***    
  theDNAElastic = new G4DNAElastic("alpha+_G4DNAElastic");
  modEI = new G4DNAIonElasticModel();
  modEI->SelectStationary(true);
  theDNAElastic->SetEmModel(modEI);
  ph->RegisterProcess(theDNAElastic, part);
      
  // *** Excitation ***
  theDNAExc = new G4DNAExcitation("alpha+_G4DNAExcitation");
  modMGE = new G4DNAMillerGreenExcitationModel();
  modMGE->SelectStationary(true);
  theDNAExc->SetEmModel(modMGE);
  ph->RegisterProcess(theDNAExc, part);
            
  // *** Ionisation ***
  theDNAIoni = new G4DNAIonisation("alpha+_G4DNAIonisation");
  modRI = new G4DNARuddIonisationModel();
  theDNAIoni->SetEmModel(new G4DNARuddIonisationModel());
  modRI->SelectStationary(true);
  ph->RegisterProcess(theDNAIoni, part);

  // *** Charge decrease ***
  theDNAChargeDecreaseProcess = 
    new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease");
  modDCD = new G4DNADingfelderChargeDecreaseModel();
  modDCD->SelectStationary(true);
  theDNAChargeDecreaseProcess->SetEmModel(modDCD);
  ph->RegisterProcess(theDNAChargeDecreaseProcess, part);

  // *** Charge increase ***
  theDNAChargeIncreaseProcess = 
    new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease");
  modDCI = new G4DNADingfelderChargeIncreaseModel();
  modDCI->SelectStationary(true);
  theDNAChargeIncreaseProcess->SetEmModel(modDCI);
  ph->RegisterProcess(theDNAChargeIncreaseProcess, part);

  // helium
  part = genericIonsManager->GetIon("helium");
      
  // *** Elastic ***    
  theDNAElastic = new G4DNAElastic("helium_G4DNAElastic");
  modEI = new G4DNAIonElasticModel();
  modEI->SelectStationary(true);
  theDNAElastic->SetEmModel(modEI);
  ph->RegisterProcess(theDNAElastic, part);
      
  // *** Excitation ***
  theDNAExc = new G4DNAExcitation("helium_G4DNAExcitation");
  modMGE = new G4DNAMillerGreenExcitationModel();
  modMGE->SelectStationary(true);
  theDNAExc->SetEmModel(modMGE);
  ph->RegisterProcess(theDNAExc, part);
            
  // *** Ionisation ***
  theDNAIoni = new G4DNAIonisation("helium_G4DNAIonisation");
  modRI = new G4DNARuddIonisationModel();
  theDNAIoni->SetEmModel(new G4DNARuddIonisationModel());
  modRI->SelectStationary(true);
  ph->RegisterProcess(theDNAIoni, part);

  // *** Charge increase ***
  theDNAChargeIncreaseProcess = 
    new G4DNAChargeIncrease("helium_G4DNAChargeIncrease");
  modDCI = new G4DNADingfelderChargeIncreaseModel();
  modDCI->SelectStationary(true);
  theDNAChargeIncreaseProcess->SetEmModel(modDCI);
  ph->RegisterProcess(theDNAChargeIncreaseProcess, part);
    
  // other ions
  part = G4GenericIon::GenericIon();

  // *** Ionisation ***
  theDNAIoni = new G4DNAIonisation("GenericIon_G4DNAIonisation");
  G4DNARuddIonisationExtendedModel* mod =
    new G4DNARuddIonisationExtendedModel();
  mod->SelectStationary(true);
  theDNAIoni->SetEmModel(mod);
  ph->RegisterProcess(theDNAIoni, part);
}

