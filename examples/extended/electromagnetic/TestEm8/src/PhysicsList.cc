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
/// \file electromagnetic/TestEm8/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// $Id: PhysicsList.cc 90770 2015-06-09 12:12:47Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsList
//
// Description: EM physics with a possibility to add PAI model
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4DecayPhysics.hh"

#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmConfigurator.hh"
#include "G4EmParameters.hh"

#include "StepMax.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList(),
  fEmPhysicsList(0),
  fDecayPhysicsList(0),
  fStepMaxProcess(0),
  fMessenger(0),
  fPAI(false)
{
  G4EmParameters::Instance()->SetVerbose(1);

  SetDefaultCutValue(1*mm);
 
  fMessenger = new PhysicsListMessenger(this);

  fStepMaxProcess = new StepMax();

  // Decay Physics is always defined
  fDecayPhysicsList = new G4DecayPhysics();

  // EM physics
  fEmName = G4String("emstandard");
  fEmPhysicsList = new G4EmStandardPhysics(1);

  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
  delete fDecayPhysicsList;
  delete fEmPhysicsList;
  for(size_t i=0; i<fHadronPhys.size(); ++i) { delete fHadronPhys[i]; }
  delete fStepMaxProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  fDecayPhysicsList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  fEmPhysicsList->ConstructProcess();
  if(fPAI) { AddPAIModel(fEmName); }
  fDecayPhysicsList->ConstructProcess();
  for(size_t i=0; i<fHadronPhys.size(); ++i) { 
    fHadronPhys[i]->ConstructProcess(); 
  }
  AddStepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) {
    return;

  } else if (name == "emstandard_opt1") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1();

  } else if (name == "emstandard_opt2") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "emstandard_opt3") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();

  } else if (name == "emstandard_opt4") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4();

  } else if (name == "emlivermore") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();

  } else if (name == "empenelope") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();

  } else if (name == "emlowenergy") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLowEPPhysics();

  } else if (name == "pai") {

    fEmName = name;
    fPAI = true;

  } else if (name == "pai_photon") { 

    fEmName = name;
    fPAI = true;

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process

  theParticleIterator->reset();
  while ((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (fStepMaxProcess->IsApplicable(*particle))
    {
      pmanager ->AddDiscreteProcess(fStepMaxProcess);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100.*eV,1e5);
  if ( verboseLevel > 0 ) { DumpCutValuesTable(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPAIModel(const G4String& modname)
{
  theParticleIterator->reset();
  while ((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String partname = particle->GetParticleName();
    if(partname == "e-" || partname == "e+") {
      NewPAIModel(particle, modname, "eIoni");

    } else if(partname == "mu-" || partname == "mu+") {
      NewPAIModel(particle, modname, "muIoni");

    } else if(partname == "proton" ||
              partname == "pi+" ||
              partname == "pi-"   
              ) {
      NewPAIModel(particle, modname, "hIoni");
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::NewPAIModel(const G4ParticleDefinition* part, 
                              const G4String& modname,
                              const G4String& procname)
{
  G4EmConfigurator* config = G4LossTableManager::Instance()->EmConfigurator();
  G4String partname = part->GetParticleName();
  if(modname == "pai") {
    G4PAIModel* pai = new G4PAIModel(part,"PAIModel");
    config->SetExtraEmModel(partname,procname,pai,"GasDetector",
                            0.0,100.*TeV,pai);
  } else if(modname == "pai_photon") {
    G4PAIPhotModel* pai = new G4PAIPhotModel(part,"PAIPhotModel");
    config->SetExtraEmModel(partname,procname,pai,"GasDetector",
                            0.0,100.*TeV,pai);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

