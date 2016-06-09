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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//
//
// Hadrontherapy Users are recommended to use the prepared macro files in order to activate the 
// most appropriate physics for hadrontherapy applications.
// As one can easily see the physics we suggest is contained in the
// 'QGSP_BIC_EMY' list adding to the decay processes (activated as default).
//
//
//    ******                     SUGGESTED PHYSICS                 *********
//
//    AT MOMENT, IF ACCURATE RESULTS ARE NEEDED, WE STRONGLY RECOMMEND: 
//    1. The use of the macro 'hadron_therapy.mac', that is connected with the HadrontherapyPhysicsList.cc file.
//    2. the QGSP_BIC_EMY Reference Physics Lists (define the PHYSLIST eviroment variable):
//       export PHYSLIST=QGSP_BIC_EMY
//       User must considered that, in this second case, radioactive processes are not activated 
 
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPhysicsListMessenger.hh"
#include "HadrontherapyStepMax.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"

// Local physic directly implemented in the Hadronthrapy directory
#include "LocalIonIonInelasticPhysic.hh"             // Physic dedicated to the ion-ion inelastic processes

// Physic lists (contained inside the Geant4 source code, in the 'physicslists folder')
#include "HadronPhysicsQGSP_BIC.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4EmProcessOptions.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyPhysicsList::HadrontherapyPhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  helIsRegistered  = false;
  bicIsRegistered  = false;
  biciIsRegistered = false;
  locIonIonInelasticIsRegistered = false;
  radioactiveDecayIsRegistered = false;

  stepMaxProcess  = 0;

  pMessenger = new HadrontherapyPhysicsListMessenger(this);

  SetVerboseLevel(1);

  // ******     Definition of some defaults for the physics     *****
  // ******     in case no physics is called by the macro file  *****
  // EM physics
  emPhysicsList = new G4EmStandardPhysics_option3(1);
  emName = G4String("emstandard_opt3");

  // Decay physics and all particles
  decPhysicsList = new G4DecayPhysics();
  raddecayList = new G4RadioactiveDecayPhysics();
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyPhysicsList::~HadrontherapyPhysicsList()
{
  delete pMessenger;
  delete emPhysicsList;
  delete decPhysicsList;
  delete raddecayList;

  for(size_t i=0; i<hadronPhys.size(); i++) {delete hadronPhys[i];}
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::ConstructParticle()
{
  decPhysicsList->ConstructParticle();
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::ConstructProcess()
{
  // transportation
  AddTransportation();

  // electromagnetic physics list
  emPhysicsList->ConstructProcess();
  em_config.AddModels();

  // decay physics list
  decPhysicsList->ConstructProcess();
  raddecayList->ConstructProcess();

  // hadronic physics lists
  for(size_t i=0; i<hadronPhys.size(); i++) {
    hadronPhys[i] -> ConstructProcess();
  }

  // step limitation (as a full process)
  //
  AddStepMax();
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::AddPhysicsList(const G4String& name)
{

  if (verboseLevel>1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }
  if (name == emName) return;

  /////////////////////////////////////////////////////////////////////////////
  //   ELECTROMAGNETIC MODELS
  /////////////////////////////////////////////////////////////////////////////
  if (name == "standard_opt3") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option3();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;

  } else if (name == "LowE_Livermore") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmLivermorePhysics();
    G4RunManager::GetRunManager()-> PhysicsHasBeenModified();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmLivermorePhysics" << G4endl;

  } else if (name == "LowE_Penelope") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmPenelopePhysics();
    G4RunManager::GetRunManager()-> PhysicsHasBeenModified();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmPenelopePhysics" << G4endl;
    
  } else if (name == "local_ion_ion_inelastic") {
    hadronPhys.push_back(new LocalIonIonInelasticPhysic());
    locIonIonInelasticIsRegistered = true;

  } else if (name == "QGSP_BIC_EMY") {
    AddPhysicsList("emstandard_opt3");
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC());
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4HadronElasticPhysics());
    hadronPhys.push_back( new G4StoppingPhysics());
    hadronPhys.push_back( new G4IonBinaryCascadePhysics());
    hadronPhys.push_back( new G4NeutronTrackingCut());
    
  } else { 
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
	   << " is not defined"
	   << G4endl;
  }
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::AddStepMax()
{
  // Step limitation seen as a process
  stepMaxProcess = new HadrontherapyStepMax();

  theParticleIterator->reset();
  while ((*theParticleIterator)()){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (stepMaxProcess->IsApplicable(*particle) && pmanager)
      {
	pmanager ->AddDiscreteProcess(stepMaxProcess);
      }
  }
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::SetCuts()
{

  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  // Set cuts for detector
  SetDetectorCut(defaultCutValue); 
  if (verboseLevel>0) DumpCutValuesTable();
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

void HadrontherapyPhysicsList::SetDetectorCut(G4double cut)
{
  G4String regionName = "DetectorLog";
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(regionName);

  G4ProductionCuts* cuts = new G4ProductionCuts ;
  cuts -> SetProductionCut(cut,G4ProductionCuts::GetIndex("gamma"));
  cuts -> SetProductionCut(cut,G4ProductionCuts::GetIndex("e-"));
  cuts -> SetProductionCut(cut,G4ProductionCuts::GetIndex("e+"));
  region -> SetProductionCuts(cuts);
}

