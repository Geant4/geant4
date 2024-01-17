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
// Code developed by
// Silvia Pozzi (1), silvia.pozzi@iss.it
// Barbara Caccia (1), barbara.caccia@iss.it
// Carlo Mancini Terracciano (2), carlo.mancini.terracciano@roma1.infn.it
// (1) Istituto Superiore di Sanita' and INFN Roma, Italy
// (2) Univ. La Sapienza and INFN Roma, Italy

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"
#include "StepMax.hh"

#include <G4EmStandardPhysics.hh>
#include <G4eBremsstrahlung.hh>

#include <G4DecayPhysics.hh> 
#include <G4ProductionCutsTable.hh>
#include <G4SystemOfUnits.hh>

#include "G4SystemOfUnits.hh"
#include "G4VPhysicsConstructor.hh"

// Physics lists
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4Decay.hh"

#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"

/////////////////////////////////////////////////////////////////////////////
PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
	defaultCutValue = 1.*mm;
	helIsRegistered  = false;
	bicIsRegistered  = false;
	biciIsRegistered = false;
	locIonIonInelasticIsRegistered = false;

	stepMaxProcess = nullptr;

	pMessenger = new PhysicsListMessenger(this);

	SetVerboseLevel(1);

	// EM physics
	emPhysicsList = new G4EmStandardPhysics_option3(1);
	emName = G4String("emstandard_opt3");

	//  emPhysicsList = new G4EmLivermorePhysics();
	//  emName = G4String("LowE_Livermore");

	// Decay physics and all particles
	decPhysicsList = new G4DecayPhysics();
}

void PhysicsList::SetCuts()
{
	G4VUserPhysicsList::SetVerboseLevel(2);
	G4VUserPhysicsList::SetCuts();
}

/////////////////////////////////////////////////////////////////////////////
PhysicsList::~PhysicsList()
{
	delete pMessenger;
	delete emPhysicsList;
	delete decPhysicsList;
	for(size_t i=0; i<hadronPhys.size(); i++) {delete hadronPhys[i];}
}

void PhysicsList::ConstructParticle()
{
	decPhysicsList->ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
	// transportation
	AddTransportation();

	// electromagnetic physics list
	emPhysicsList->ConstructProcess();

	// decay physics list
	decPhysicsList->ConstructProcess();

	// hadronic physics lists
	for(size_t i=0; i<hadronPhys.size(); i++)
	{
		hadronPhys[i]->ConstructProcess();
	}

	// step limitation (as a full process)
	AddStepMax();
}

void PhysicsList::AddPhysicsList(const G4String& name)
{
	if (verboseLevel>1)
	{
		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
	}
	if (name == emName) return;

	//   Electromagnetic Models
	if (name == "standard_opt3")
	{
		emName = name;
		delete emPhysicsList;
		emPhysicsList = new G4EmStandardPhysics_option3();
		G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: "
				<< "G4EmStandardPhysics_option3" << G4endl;
	}
	else if (name == "standard_opt4")
	{
		emName = name;
		delete emPhysicsList;
		emPhysicsList = new G4EmStandardPhysics_option4();
		G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: "
				<< "G4EmStandardPhysics_option4" << G4endl;
	}
	else if (name == "LowE_Livermore")
	{
		emName = name;
		delete emPhysicsList;
		emPhysicsList = new G4EmLivermorePhysics();
		G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: "
				<< "G4EmLivermorePhysics" << G4endl;
	}
	else if (name == "LowE_Penelope")
	{
		emName = name;
		delete emPhysicsList;
		emPhysicsList = new G4EmPenelopePhysics();
		G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: "
				<< "G4EmLivermorePhysics" << G4endl;

	//   Hadronic Models
	}
	else if (name == "elastic" && !helIsRegistered)
	{

		hadronPhys.push_back( new G4HadronElasticPhysics());
		helIsRegistered = true;
		G4cout << "THE FOLLOWING HADRONIC ELASTIC PHYSICS LIST HAS BEEN ACTIVATED: "
				<< "G4HadronElasticPhysics()" << G4endl;
	}
	else if (name == "binary" && !bicIsRegistered)
	{
		hadronPhys.push_back(new G4HadronInelasticQBBC());
		bicIsRegistered = true;
		G4cout << "THE FOLLOWING HADRONIC INELASTIC PHYSICS LIST HAS BEEN ACTIVATED: "
				<< "G4HadronInelasticQBBC()" << G4endl;
	}
	else if (name == "binary_ion" && !biciIsRegistered)
	{
		hadronPhys.push_back(new G4IonBinaryCascadePhysics());
		biciIsRegistered = true;
	}
	else
	{
		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
			   << " is not defined"
			   << G4endl;
	}
}

void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process
	stepMaxProcess = new StepMax();

	auto particleIterator=GetParticleIterator();
	particleIterator->reset();
	while ((*particleIterator)())
	{
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (stepMaxProcess->IsApplicable(*particle) && pmanager)
    {
    	pmanager ->AddDiscreteProcess(stepMaxProcess);
    }
  }
}
