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
// HadrontherapyPhysicsList.cc
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy

// This class provides all the physic models that can be activated inside Hadrontherapy;
// Each model can be setted via macro commands;
// Inside Hadrontherapy the models can be activate with three different complementar methods:
//
//    For Hadrontherapy we suggest the use of:
//
//    /Physic/addPhysics/emstandard_option3 (electromagnetic model)
//    /Physic/addPhysics/QElastic (hadronic elastic model)
//    /Physic/addPhysics/binary (hadronic inelastic models for proton and neutrons)
//    /Physic/addPhysics/binary_ion (hadronic inelastic models for ions)
//
// or alternatevely write only the following command in the macro file to activate
// a Reference physics list containing all the processes:
//
//    /Physics/addPhysics QGSP_BIC_EMY 
//
//    Example of the use of physics lists can be found in the macro files included in the
//    'macro' folder .
//
// 3. Use of a *local* physics. In this case the models are implemented in local files
//    contained in the Hadrontherapy folder. The use of local physic is recommended
//    to more expert Users.
//    We provide as local, only the LocalStandardICRU73EmPhysic.cc (an Elecromagnetic
//    implementation containing the new ICRU73 data table for ions stopping powers)
//    and the LocalIonIonInelasticPhysic.cc (physic list to use for the ion-ion interaction
//    case)
//    The *local* physics can be activated with the same /Physic/addPhysics <nameOfPhysic> command;
//
//    AT MOMENT, IF ACCURATE RESULTS ARE NEDED, WE STRONGLY RECOMMEND THE USE OF THE MACROS:
//    proton_therapy.mac: for proton beams;
//    ion_therapy.mac   : for ion beams;
 
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
#include "LocalINCLIonIonInelasticPhysic.hh"         // Physic dedicated to the ion-ion inelastic processes using INCL/ABLA

// Physic lists (contained inside the Geant4 source code, in the 'physicslists folder')
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"

#include "G4QStoppingPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4HadronQElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4Decay.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "HadronPhysicsQGSP_BIC.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4EmProcessOptions.hh"

#include "G4RadioactiveDecayPhysics.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyPhysicsList::HadrontherapyPhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  helIsRegisted  = false;
  bicIsRegisted  = false;
  biciIsRegisted = false;
  locIonIonInelasticIsRegistered = false;
  radioactiveDecayIsRegisted = false;

  stepMaxProcess  = 0;

  pMessenger = new HadrontherapyPhysicsListMessenger(this);

  SetVerboseLevel(1);

  // EM physics
  emPhysicsList = new G4EmStandardPhysics_option3(1);
  emName = G4String("emstandard_opt3");

  // Deacy physics and all particles
  decPhysicsList = new G4DecayPhysics();
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyPhysicsList::~HadrontherapyPhysicsList()
{
  delete pMessenger;
  delete emPhysicsList;
  delete decPhysicsList;
  for(size_t i=0; i<hadronPhys.size(); i++) {delete hadronPhys[i];}
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::AddPackage(const G4String& name)
{
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys =factory.GetReferencePhysList(name);
  G4int i=0;
  const G4VPhysicsConstructor* elem= phys->GetPhysics(i);
  G4VPhysicsConstructor* tmp = const_cast<G4VPhysicsConstructor*> (elem);
  while (elem !=0)
    {
      RegisterPhysics(tmp);
      elem= phys->GetPhysics(++i) ;
      tmp = const_cast<G4VPhysicsConstructor*> (elem);
    }
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

  // hadronic physics lists
  for(size_t i=0; i<hadronPhys.size(); i++) {
    hadronPhys[i]->ConstructProcess();
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
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmLivermorePhysics" << G4endl;

    /////////////////////////////////////////////////////////////////////////////
    //   HADRONIC MODELS
    /////////////////////////////////////////////////////////////////////////////
  } else if (name == "elastic" && !helIsRegisted) {
    G4cout << "THE FOLLOWING HADRONIC ELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4HadronElasticPhysics()" << G4endl;
    hadronPhys.push_back( new G4HadronElasticPhysics());
    helIsRegisted = true;

  } else if (name == "DElastic" && !helIsRegisted) {
    hadronPhys.push_back( new G4HadronDElasticPhysics());
    helIsRegisted = true;

  } else if (name == "HElastic" && !helIsRegisted) {
    hadronPhys.push_back( new G4HadronHElasticPhysics());
    helIsRegisted = true;

  } else if (name == "QElastic" && !helIsRegisted) {
    hadronPhys.push_back( new G4HadronQElasticPhysics());
    helIsRegisted = true;

  } else if (name == "binary" && !bicIsRegisted) {
    hadronPhys.push_back(new G4HadronInelasticQBBC());
    bicIsRegisted = true;
    G4cout << "THE FOLLOWING HADRONIC INELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4HadronInelasticQBBC()" << G4endl;

  } else if (name == "binary_ion" && !biciIsRegisted) {
    hadronPhys.push_back(new G4IonBinaryCascadePhysics());
    biciIsRegisted = true;

  } else if (name == "local_ion_ion_inelastic" && !locIonIonInelasticIsRegistered) {
    hadronPhys.push_back(new LocalIonIonInelasticPhysic());
    locIonIonInelasticIsRegistered = true;

  } else if (name == "local_incl_ion_ion_inelastic" && !locIonIonInelasticIsRegistered) {
    hadronPhys.push_back(new LocalINCLIonIonInelasticPhysic());
    locIonIonInelasticIsRegistered = true;
      
  } else if (name == "radioactive_decay" && !radioactiveDecayIsRegisted ) {
    hadronPhys.push_back(new G4RadioactiveDecayPhysics());
    radioactiveDecayIsRegisted = true;
      
  } else if (name == "QGSP_BIC_EMY") {
    AddPhysicsList("emstandard_opt3");
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4HadronElasticPhysics());
    hadronPhys.push_back( new G4QStoppingPhysics());
    hadronPhys.push_back( new G4IonBinaryCascadePhysics());
    hadronPhys.push_back( new G4NeutronTrackingCut());
    hadronPhys.push_back( new HadronPhysicsQGSP_BIC());
      
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

