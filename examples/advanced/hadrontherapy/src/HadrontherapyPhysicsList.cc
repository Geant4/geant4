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
// Physics models in Hadrontherapy, following the Geant4 organisation, can be definided using three different approaches:
// 1. Activating one of the 'Reference Physics Lists' that are already prepared by
//    the Geant4 Collaboration and are contained in the $G4INSTALL/source/physics_lists/lists folder
//    The 'Reference Physics Lists' can be activated setting a specific enviroment variable to the name
//    of the physics. For example if the QGSP_BIC Reference Physics Lists must be activated the User 
//    must set export PHYSLIST=QGSP_BIC (or setenv PHYSLIST QGSP_BIC).
//    A 'Reference Physics Lists' contains all the physics process necessary to a particle transport
//    If the User set the PHYSLIST variable Hadrontherapy will start with the defaultMacroWithReferencePhysicsList.mac
//    macro. See this macro file for more details
//
// 2. Activating the 'Builders' already prepared by
//    the Geant4 Collaboration and contained in the $G4INSTALL/source/physics_lists/builder folder.
//    Each builder is specific of a given model. There are builders for the electromagnetic processes, for the
//    hadronic one, etc.
//    If the PHYSLIST variable is not defined Hadrontherapy starts with the defaultMacro.mac where the single builders
//    are activated for the various processes of interest.
//    Each builder is activated with the /Physics/addPhysics <nome builder> command
//
// 3. Defining a specific 'local' physics list. In Hadrontherapy two loca physics list are defined 
//    (LocalINCLIonIonInelasticPhysic.cc and LocalIonIonInelasticPhysic.cc) where models for the 
//    ion-ion interactions are activated.
//    'Local' physics can be activated using the /Physics/addPhysics <nome builder> command
//    (see the defaultMacro.mac to se an example).
//
//    ******       SUGGESTED PHYSICS       *********
//
//    AT MOMENT, IF ACCURATE RESULTS ARE NEDED, WE STRONGLY RECOMMEND: 
//    1. The use of the macro 'hadron_therapy.mac', or
//    2. the QGSP_BIC_EMY Reference Physics Lists (define the PHYSLIST eviroment variable):
//       export PHYSLIST=QGSP_BIC_EMY
 
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
#include "G4DecayPhysics.hh"
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

  helIsRegistered  = false;
  bicIsRegistered  = false;
  biciIsRegistered = false;
  locIonIonInelasticIsRegistered = false;
  radioactiveDecayIsRegistered = false;

  stepMaxProcess  = 0;

  pMessenger = new HadrontherapyPhysicsListMessenger(this);

  SetVerboseLevel(1);

  // EM physics
  emPhysicsList = new G4EmStandardPhysics_option3(1);
  emName = G4String("emstandard_opt3");

  // Decay physics and all particles
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
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmLivermorePhysics" << G4endl;

	/////////////////////////////////////////////////////////////////////////////
	//   HADRONIC MODELS
	/////////////////////////////////////////////////////////////////////////////
    } else if (name == "Elastic")
    {
	if(!helIsRegistered) 
	{
	    G4cout << "THE FOLLOWING HADRONIC ELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4HadronElasticPhysics()" << G4endl;
	    hadronPhys.push_back( new G4HadronElasticPhysics());
	    helIsRegistered = true;
	}
	else  G4cout << "AN ELASTIC PHYSICS HAS BEEN ALREADY ACTIVATED!" << G4endl;
    }
    else if (name == "DElastic")
    {
	if(!helIsRegistered) 
	{
	    hadronPhys.push_back( new G4HadronDElasticPhysics());
	    helIsRegistered = true;
	}
	else  G4cout << "AN ELASTIC PHYSICS HAS BEEN ALREADY ACTIVATED!" << G4endl;

    }
    else if (name == "HElastic")
    {
	if(!helIsRegistered) 
	{
	    hadronPhys.push_back( new G4HadronHElasticPhysics());
	    helIsRegistered = true;
	}
	else  G4cout << "AN ELASTIC PHYSICS HAS BEEN ALREADY ACTIVATED!" << G4endl;

    }
    else if (name == "QElastic")
    {
	if(!helIsRegistered) 
	{
	    hadronPhys.push_back( new G4HadronQElasticPhysics());
	    helIsRegistered = true;
	}
	else  G4cout << "AN ELASTIC PHYSICS HAS BEEN ALREADY ACTIVATED!" << G4endl;

    }
    else if (name == "Em_extra_physics")
    {
	    hadronPhys.push_back( new G4EmExtraPhysics());
    }
    else if (name == "Stopping_physics")
    {
	    hadronPhys.push_back( new G4QStoppingPhysics());
    }
    else if (name == "Neutron_tracking_cut")
    {
	    hadronPhys.push_back( new G4NeutronTrackingCut());
    }
    else if (name == "Hadron_QGSP_BIC")
    {
	    hadronPhys.push_back( new HadronPhysicsQGSP_BIC());
	   // helIsRegistered = true;
    }
    else if (name == "Hadron_QBBC") 
    {
	    hadronPhys.push_back(new G4HadronInelasticQBBC());
	    //bicIsRegistered = true;
	    G4cout << "THE FOLLOWING HADRONIC INELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4HadronInelasticQBBC()" << G4endl;
    } 

    else if (name == "binary") 
    {
    hadronPhys.push_back(new G4HadronInelasticQBBC());
    //bicIsRegisted = true;
    G4cout << "THE FOLLOWING HADRONIC INELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4HadronInelasticQBBC()" << G4endl;
    }
    
    else if (name == "binary_ion")
    {
	    hadronPhys.push_back(new G4IonBinaryCascadePhysics());
	    //biciIsRegistered = true;
    }
    else if (name == "local_ion_ion_inelastic")
    {
	    hadronPhys.push_back(new LocalIonIonInelasticPhysic());
	    locIonIonInelasticIsRegistered = true;
    }
    else if (name == "local_incl_ion_ion_inelastic")
    {
	    hadronPhys.push_back(new LocalINCLIonIonInelasticPhysic());
	    locIonIonInelasticIsRegistered = true;
    } 
    else if (name == "decay")
    {
	hadronPhys.push_back(new G4DecayPhysics());
	//radioactiveDecayIsRegistered = true;
    }
    else if (name == "radioactive_decay" && !radioactiveDecayIsRegistered )
    {
	hadronPhys.push_back(new G4RadioactiveDecayPhysics());
	radioactiveDecayIsRegistered = true;

	// The following is the construction of the QGSP_BIC_EMY Reference physics list 
	// reconstructed here like a builder: it should be identical to the
	// one contained inside the $G4INSTALL/physics_lists/lists folder
    }
    else if (name == "QGSP_BIC_EMY")
    {
	AddPhysicsList("emstandard_opt3");
	hadronPhys.push_back( new G4EmExtraPhysics());
	hadronPhys.push_back( new G4HadronElasticPhysics());
	hadronPhys.push_back( new G4QStoppingPhysics());
	hadronPhys.push_back( new G4IonBinaryCascadePhysics());
	hadronPhys.push_back( new G4NeutronTrackingCut());
	hadronPhys.push_back( new HadronPhysicsQGSP_BIC());
	hadronPhys.push_back( new G4DecayPhysics());

    }
    else {

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

