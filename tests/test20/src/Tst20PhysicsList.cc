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
// $Id: Tst20PhysicsList.cc,v 1.14 2010-04-01 08:17:10 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------

#include "Tst20PhysicsList.hh"
#include "G4DNAGenericIonsManager.hh"


#include "Tst20PhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"

#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"

// e+
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"

#include "G4Decay.hh"

#include "G4FinalStateProduct.hh"
#include "G4DNAProcess.hh"

#include "G4CrossSectionExcitationEmfietzoglou.hh"
#include "G4FinalStateExcitationEmfietzoglou.hh"

#include "G4CrossSectionElasticScreenedRutherfordLE.hh"
#include "G4CrossSectionElasticScreenedRutherfordHE.hh"
#include "G4FinalStateElasticScreenedRutherford.hh"
#include "G4FinalStateElasticBrennerZaider.hh"

#include "G4CrossSectionElasticChampion.hh"
#include "G4FinalStateElasticChampion.hh"

#include "G4CrossSectionExcitationBorn.hh"
#include "G4FinalStateExcitationBorn.hh"

#include "G4CrossSectionIonisationBorn.hh"
#include "G4FinalStateIonisationBorn.hh"

#include "G4CrossSectionIonisationRudd.hh"
#include "G4FinalStateIonisationRudd.hh"

#include "G4CrossSectionExcitationMillerGreen.hh"
#include "G4FinalStateExcitationMillerGreen.hh"

#include "G4CrossSectionChargeDecrease.hh"
#include "G4FinalStateChargeDecrease.hh"

#include "G4CrossSectionChargeIncrease.hh"
#include "G4FinalStateChargeIncrease.hh"

#include "G4CrossSectionKill.hh"
#include "G4FinalStateKill.hh"

// Processes

typedef G4DNAProcess<G4CrossSectionElasticScreenedRutherfordHE,G4FinalStateElasticScreenedRutherford> ElasticScreenedRutherfordHE;
typedef G4DNAProcess<G4CrossSectionElasticScreenedRutherfordLE,G4FinalStateElasticBrennerZaider> ElasticScreenedRutherfordLE;
typedef G4DNAProcess<G4CrossSectionElasticChampion,G4FinalStateElasticChampion> ElasticChampion;
typedef G4DNAProcess<G4CrossSectionExcitationEmfietzoglou,G4FinalStateExcitationEmfietzoglou> ExcitationEmfietzoglou;
typedef G4DNAProcess<G4CrossSectionExcitationBorn,G4FinalStateExcitationBorn> ExcitationBorn;
typedef G4DNAProcess<G4CrossSectionIonisationBorn,G4FinalStateIonisationBorn> IonisationBorn;
typedef G4DNAProcess<G4CrossSectionIonisationRudd,G4FinalStateIonisationRudd> IonisationRudd;
typedef G4DNAProcess<G4CrossSectionExcitationMillerGreen,G4FinalStateExcitationMillerGreen> ExcitationMillerGreen;
typedef G4DNAProcess<G4CrossSectionChargeDecrease,G4FinalStateChargeDecrease> ChargeDecrease;
typedef G4DNAProcess<G4CrossSectionChargeIncrease,G4FinalStateChargeIncrease> ChargeIncrease;
typedef G4DNAProcess<G4CrossSectionKill,G4FinalStateKill> KillBelowThreshold;


Tst20PhysicsList::Tst20PhysicsList(): G4VUserPhysicsList()
{
  defaultCutValue = 1000. * nanometer;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  cutForProton    = defaultCutValue;
  // SetVerboseLevel(1);
  physicsListMessenger = new Tst20PhysicsListMessenger(this);
}


Tst20PhysicsList::~Tst20PhysicsList()
{
  delete physicsListMessenger;
}


void Tst20PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
}


void Tst20PhysicsList::ConstructBosons()
{
  
  // gamma
  G4Gamma::GammaDefinition();
  
}

void Tst20PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}


void Tst20PhysicsList::ConstructBarions()
{
  //  baryons
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  // Geant4-DNA

  G4DNAGenericIonsManager* genericIonsManager;
  genericIonsManager = G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
}


void Tst20PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}


void Tst20PhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* processManager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      
      if (particleName == "gamma") 
	{
	  processManager->AddDiscreteProcess(new G4LowEnergyCompton);
	  processManager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
	  processManager->AddDiscreteProcess(new G4LowEnergyPhotoElectric);	
	  processManager->AddDiscreteProcess(new G4LowEnergyRayleigh);	  
	} 
      else if (particleName == "e-") 
	{
	  // DNA
	  processManager->AddDiscreteProcess(new ExcitationEmfietzoglou("ExcitationEmfietzoglou"));
/*
          processManager->AddDiscreteProcess(new ElasticScreenedRutherfordLE("ElasticScreenedRutherfordLE"));
          processManager->AddDiscreteProcess(new ElasticScreenedRutherfordHE("ElasticScreenedRutherfordHE"));
*/
          processManager->AddDiscreteProcess(new ElasticChampion("ElasticChampion"));
	  processManager->AddDiscreteProcess(new IonisationBorn("IonisationBorn"));

	}
      else if (particleName == "e+") 
	{
	  processManager->AddProcess(new G4eMultipleScattering,-1, 1,1);
	  processManager->AddProcess(new G4eIonisation,      -1, 2,2);
	  processManager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);
	  processManager->AddProcess(new G4eplusAnnihilation,  0,-1,4);	
	} 
      else if (particleName == "proton") 
	{
	  processManager->AddDiscreteProcess(new ExcitationMillerGreen);
	  processManager->AddDiscreteProcess(new ExcitationBorn);
	  processManager->AddDiscreteProcess(new IonisationRudd);
	  processManager->AddDiscreteProcess(new IonisationBorn);
	  processManager->AddDiscreteProcess(new ChargeDecrease);
	} 
      else if (particleName == "hydrogen") 
	{
	  processManager->AddDiscreteProcess(new IonisationRudd);
	  processManager->AddDiscreteProcess(new ChargeIncrease);
	} 
      else if (particleName == "alpha") 
	{
	  processManager->AddDiscreteProcess(new ExcitationMillerGreen);
	  processManager->AddDiscreteProcess(new IonisationRudd);
	  processManager->AddDiscreteProcess(new ChargeDecrease);   
	} 
      else if (particleName == "alpha+") 
	{
	  processManager->AddDiscreteProcess(new ExcitationMillerGreen);
	  processManager->AddDiscreteProcess(new IonisationRudd);
	  processManager->AddDiscreteProcess(new ChargeDecrease);
	  processManager->AddDiscreteProcess(new ChargeIncrease);    
	} 
      else if (particleName == "helium") 
	{
	  processManager->AddDiscreteProcess(new ExcitationMillerGreen);
	  processManager->AddDiscreteProcess(new IonisationRudd);
	  processManager->AddDiscreteProcess(new ChargeIncrease);
	}
    }
}

void Tst20PhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* processManager = particle->GetProcessManager();
      if (theDecayProcess->IsApplicable(*particle)) 
	{ 
	  processManager ->AddProcess(theDecayProcess);
	  // set ordering for PostStepDoIt and AtRestDoIt
	  processManager ->SetProcessOrdering(theDecayProcess, idxPostStep);
	  processManager ->SetProcessOrdering(theDecayProcess, idxAtRest);
	}
    }
}


void Tst20PhysicsList::SetCuts()
{
  if (verboseLevel >0)
  {
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 
  SetCutValue(cutForProton, "proton");
  SetCutValue(cutForProton, "anti_proton");
  
  if (verboseLevel>0) DumpCutValuesTable();
}



