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
// -------------------------------------------------------------------
// $Id: PhysicsList.cc,v 1.3 2010-09-09 10:16:15 sincerti Exp $
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "PhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsList::PhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1*nanometer;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  cutForProton    = defaultCutValue;
  
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsList::~PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructBosons()
{ 

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//DNA
#include "G4DNAGenericIonsManager.hh"
//ENDDNA

void PhysicsList::ConstructBarions()
{
  //  baryons
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  // Geant4 DNA new particles
  G4DNAGenericIonsManager * genericIonsManager;
  genericIonsManager=G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// MIGRATED Geant4-DNA MODELS

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
  {

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

// MIGRATED Geant4-DNA

    if (particleName == "e-") {

      G4DNAElastic* theDNAElasticProcess = new G4DNAElastic();
      //theDNAElasticProcess->SetModel(new G4DNAChampionElasticModel());
      theDNAElasticProcess->SetModel(new G4DNAScreenedRutherfordElasticModel());
      pmanager->AddDiscreteProcess(theDNAElasticProcess);

      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAVibExcitation());
      pmanager->AddDiscreteProcess(new G4DNAAttachment()); 

     	    
    } else if ( particleName == "proton" ) {

      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeDecrease());

    } else if ( particleName == "hydrogen" ) {

      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeIncrease());

    } else if ( particleName == "alpha" ) {

      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeDecrease());

    } else if ( particleName == "alpha+" ) {

      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeDecrease());
      pmanager->AddDiscreteProcess(new G4DNAChargeIncrease());
	    
    } else if ( particleName == "helium" ) {

      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeIncrease());

    }

  } // Loop on particles


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructGeneral()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetCuts()
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
