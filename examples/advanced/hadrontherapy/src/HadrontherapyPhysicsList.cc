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
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F.Romano(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// 
// * cirrone@lns.infn.it
//
// See more at: http://workgroup.lngs.infn.it/geant4lns/
// ----------------------------------------------------------------------------

// This class provides all the physic models that can be activated inside Hadrontherapy;
// Each model can be setted via macro commands;
// Inside Hadrontherapy the models can be activate with three different complementar methods:
//
// 1. Use of the *Packages*. 
//    Packages (that are contained inside the
//    Geant4 distribution at $G4INSTALL/source/physics_lists/lists) provide a full set 
//    of models (both electromagnetic and hadronic).
//    The User can use this method simply add the line /physic/addPackage <nameOfPackage> 
//    in his/her macro file. No other action is required.
//    For Hadrontherapy applications we suggest the use of the QGSP_BIC package 
//    for proton beams. The same can be used 
//    also for ligth ion beam.
//    Example of use of package can be found in the packageQGSP_BIC.mac file.
//
// 2. Use of the *Physic Lists*.
//    Physic lists are also already ready to use inside the Geant4 distribution 
//    ($G4INSTALL/source/physics_lists/builders). To use them the simple
//    /physic/addPhysics <nameOfPhysicList> command must be used in the macro.
//    In Hadrontherapy we provide physics list to activate Electromagnetic,
//    Hadronic elastic and Hadronic inelastic models.
//
//    For Hadrontherapy we suggest the use of:
//    
//    /physic/addPhysic/emstandard_option3 (electromagnetic model)
//    /physic/addPhysic/QElastic (hadronic elastic model)
//    /physic/addPhysic/binary (hadronic inelastic models for proton and neutrons)
//    /physic/addPhysic/binary_ion (hadronic inelastic models for ions)
//
//    Example of the use of physics lists can be found in the macro files included in the
//    'macro' folder .   
//
// 3. Use of *local* physics. In this case the models are implemented in local files
//    contained in the Hadrontherapy folder. The use of local physic is recommended 
//    to more expert Users. Actually we provide a local physic wih the implementation of the 
//    Low Energy Livermore electromagnetic models, Low Energy Penelope electromagnetic models, 
//    a local electromagnetic physic for User want the test the first implementation of 
//    ICRU73 data table in Geant4 and an class specific with the implementation of 
//    ion-ion inelastic models.
//    The *local* physic can be activated with the same /physic/addPhysic <nameOfPhysic> command;
//
//    While Packages approch must be used exclusively, Physics List and Local physics can 
//    be activated, if necessary, contemporaneously in the same simulation run.
//
//    AT MOMENT, IF ACCURATE RESULTS ARE NEDED, WE STRONGLY RECOMMEND THE USE OF THE MACROS:
//    proton_therapy.mac: use of the built-in Geant4 physics list for proton beams)
//    ion_therapy.mac   : use of mixed combination of native Geant4 physic lists 
//                        and local physic for ion-ion enelastic processes)

#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPhysicsListMessenger.hh"
#include "HadrontherapyStepMax.hh"
#include "G4PhysListFactory.hh"
//#include "G4VModularPhyisicsList.hh"
#include "G4VPhysicsConstructor.hh"

// Local physic directly implemented in the Hadronthrapy directory
#include "LocalLivermoreEmPhysic.hh"                 // Electromagnetic "Livermore" models (from the Lowenergy physic)
#include "LocalPenelopeEmPhysic.hh"                  // Electromagnetic "Penelope"  models (from the Lowenergy physic)
#include "LocalStandardICRU73EmPhysic.hh"            // This permits the use of the ICRU73 tables for stopping powers of ions
#include "LocalIonIonInelasticPhysic.hh"             // Physic dedicated to the ion-ion inelastic processes

// Physic lists (contained inside the Geant4 distribution)
#include "G4EmStandardPhysics_option3.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4HadronQElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4Decay.hh"

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
  defaultCutValue = 0.01 *mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  helIsRegisted  = false;
  bicIsRegisted  = false;
  biciIsRegisted = false;

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
  //
  AddTransportation();
  
  // electromagnetic physics list
  //
  emPhysicsList->ConstructProcess();
  em_config.AddModels();

  // decay physics list
  //
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

  if (name == "standard_opt3") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option3();    
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;
    
  } else if (name == "local_standardICRU73") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new LocalStandardICRU73EmPhysic(name);
    em_config.SetExtraEmModel("GenericIon","ionIoni",
			      new G4IonParametrisedLossModel(),
			      "",0.0, 100.0*TeV,
			      new G4IonFluctuations());
    G4cout << "standardICRU73" << G4endl;
    
  } else if (name == "local_livermore") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new LocalLivermoreEmPhysic();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;

  } else if (name == "local_penelope") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new LocalPenelopeEmPhysic();

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


