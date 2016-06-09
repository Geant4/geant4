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
// $Id$
//
// Author: Susanna Guatelli, susanna@uow.edu.au

#include "RemSimPhysicsList.hh"
#include "RemSimPhysicsListMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"

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
#include "HadronPhysicsQGSP_BIC_HP.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "HadronPhysicsQGSP_BIC.hh"

// The electromagnetic physics and the decay are 
// registered by default. The user has to execute the
// macro physics.mac to activate the hadronic component of the physics.
// Look vehicle1.mac, vehicle2.mac and moon.mac for example

RemSimPhysicsList::RemSimPhysicsList(): G4VModularPhysicsList()
					
{
  defaultCutValue = 0.1* mm;
  //  SetVerboseLevel(1);

  helIsRegistered  = false;
  bicIsRegistered  = false;
  bicIonIsRegistered = false;
  radioactiveDecayIsRegistered = false;

  messenger = new RemSimPhysicsListMessenger(this);

  SetVerboseLevel(1);

  // EM physics
  emPhysicsList = new G4EmStandardPhysics_option3(1);
 
  // Decay physics 
  decPhysicsList = new G4DecayPhysics();
  
  G4cout<< ">> PhysicsList:: G4EmStandardPhysics_option3 activated << "<<G4endl;
  G4cout<< ">> PhysicsList:: Decay activated << "<<G4endl;
}


RemSimPhysicsList::~RemSimPhysicsList()
{
  delete messenger;
  delete emPhysicsList;
  delete decPhysicsList;
  for(size_t i=0; i<hadronPhys.size(); i++) {delete hadronPhys[i];}
}

/////////////////////////////////////////////////////////////////////////////
void RemSimPhysicsList::ConstructParticle()
{
  decPhysicsList->ConstructParticle();
}

void RemSimPhysicsList::ConstructProcess()
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
}

void RemSimPhysicsList::AddPhysicsList(const G4String& name)
{

  if (verboseLevel>1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }


 if (name == "elastic" && !helIsRegistered) {
    G4cout << "THE FOLLOWING HADRONIC ELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4HadronElasticPhysics()" << G4endl;
    hadronPhys.push_back( new G4HadronElasticPhysics());
    helIsRegistered = true;

  } else if (name == "DElastic" && !helIsRegistered) {
    G4cout << "THE FOLLOWING HADRONIC ELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4HadronDElasticPhysics()" << G4endl;
    hadronPhys.push_back( new G4HadronDElasticPhysics());
    helIsRegistered = true;

  } else if (name == "HElastic" && !helIsRegistered) {
    G4cout << "THE FOLLOWING HADRONIC ELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4HadronHElasticPhysics()" << G4endl;
    hadronPhys.push_back( new G4HadronHElasticPhysics());
    helIsRegistered = true;

  } else if (name == "QElastic" && !helIsRegistered) {
    G4cout << "THE FOLLOWING HADRONIC ELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4HadronQElasticPhysics()" << G4endl;
    hadronPhys.push_back( new G4HadronQElasticPhysics());
    helIsRegistered = true;

  } else if (name == "binary" && !bicIsRegistered) {
    hadronPhys.push_back(new HadronPhysicsQGSP_BIC);
    bicIsRegistered = true;
    G4cout << "THE FOLLOWING HADRONIC INELASTIC PHYSICS LIST HAS BEEN ACTIVATED: QGSP_BIC" << G4endl;

  } else if (name == "binary_ion" && !bicIonIsRegistered) {
    hadronPhys.push_back(new G4IonBinaryCascadePhysics());
    bicIonIsRegistered = true;
    G4cout << "THE FOLLOWING HADRONIC INELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4IonBinaryCascadePhysics()" << G4endl;
  } else if (name == "radioactive_decay" && !radioactiveDecayIsRegistered ) {
    hadronPhys.push_back(new G4RadioactiveDecayPhysics());
    radioactiveDecayIsRegistered = true;
    G4cout << "THE FOLLOWING HADRONIC INELASTIC PHYSICS LIST HAS BEEN ACTIVATED: G4RadioactiveDecayPhysics()" << G4endl;
  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

void RemSimPhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCutsWithDefault();
  if (verboseLevel>0) DumpCutValuesTable();
}










