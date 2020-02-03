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
/// \file electromagnetic/TestEm10/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsList
//
// Description: EM physics with a possibility to add transition radiation
//
// Author:      V.Ivanchenko 16.12.2015
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4DecayPhysics.hh"
#include "TransitionRadiationPhysics.hh"

#include "G4EmParameters.hh"

#include "StepMax.hh"

#include "G4ProcessManager.hh"
//#include "G4ParticleTypes.hh"
//#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal StepMax* PhysicsList::fStepMaxProcess = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(DetectorConstruction* ptr) 
: G4VModularPhysicsList(),
  fVerbose(1)
{
  SetDefaultCutValue(1*CLHEP::mm);
  SetVerboseLevel(fVerbose);
 
  fMessenger = new PhysicsListMessenger(this);

  // Decay Physics is always defined to define all particles
  fDecayPhysicsList = new G4DecayPhysics(fVerbose);

  // EM physics
  fEmName = G4String("emstandard_opt0");
  fEmPhysicsList = new G4EmStandardPhysics(fVerbose);

  // Transition radiation physics
  fXTRPhysicsList = new TransitionRadiationPhysics(fVerbose, ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
  delete fDecayPhysicsList;
  delete fXTRPhysicsList;
  delete fEmPhysicsList;
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
  fDecayPhysicsList->ConstructProcess();
  fXTRPhysicsList->ConstructProcess(); 
  AddStepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (fVerbose > 1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) {
    return;

  } else if (name == "emstandard_opt0") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics(fVerbose);

  } else if (name == "emstandard_opt1") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1(fVerbose);

  } else if (name == "emstandard_opt2") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2(fVerbose);

  } else if (name == "emstandard_opt3") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3(fVerbose);

  } else if (name == "emstandard_opt4") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4(fVerbose);

  } else if (name == "emstandardWVI") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsWVI(fVerbose);

  } else if (name == "emstandardSS") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsSS(fVerbose);

  } else if (name == "emstandardGS") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsGS(fVerbose);

  } else if (name == "emlivermore") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics(fVerbose);

  } else if (name == "empenelope") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics(fVerbose);

  } else if (name == "emlowenergy") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLowEPPhysics(fVerbose);

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetXTRModel(const G4String& name)
{
  fXTRPhysicsList->SetXTRModel(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process
  fStepMaxProcess = new StepMax();

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)())
  {
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (fStepMaxProcess->IsApplicable(*particle))
    {
      pmanager->AddDiscreteProcess(fStepMaxProcess);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

