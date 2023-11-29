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
/// \file XraySPOPhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "XraySPOPhysicsList.hh"
#include "XraySPOPhysicsListMessenger.hh"

#include "G4EmStandardPhysics_SpacePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// particles
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Decay.hh"
#include "G4PhysicsListHelper.hh"
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XraySPOPhysicsList::XraySPOPhysicsList()
{
  fMessenger = new XraySPOPhysicsListMessenger(this);

  SetVerboseLevel(1);

  // EM physics
  fEmPhysicsList = new G4EmStandardPhysics();

  G4LossTableManager::Instance();
  SetDefaultCutValue(1.*mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XraySPOPhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ProcessManager.hh"

void XraySPOPhysicsList::ConstructProcess()
{
  AddTransportation();

  // electromagnetic physics list
  fEmPhysicsList->ConstructProcess();

  // Em options
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetBuildCSDARange(true);
  param->SetMscThetaLimit(0);

  AddDecay();

  AddRadioactiveDecay();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XraySPOPhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0)
  {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName)
  {
    return;
  }
  if (name == "SpacePhysics")
  {
    delete fEmPhysicsList;
    G4cout << "X - Defining SpacePhysics, no hadronic part" << G4endl;
    fEmName = name;
    fEmPhysicsList = new G4EmStandardPhysics_SpacePhysics();
  }
  else if (name == "emstandard_opt0")
  {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();
  }
  else if (name == "emstandard_opt1")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1();
  }
  else if (name == "emstandard_opt2")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2();
  }
  else if (name == "emstandard_opt3")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();
  }
  else if (name == "emstandard_opt4")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4();
  }
  else if (name == "emstandardSS")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsSS();
  }
  else if (name == "emstandardGS")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsGS();
  }
  else if (name == "emstandardWVI")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsWVI();
  }
  else if (name == "emlowenergy")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLowEPPhysics();
  }
  else if (name == "emlivermore")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();
  }
  else if (name == "empenelope")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();
  }
  else
  {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XraySPOPhysicsList::AddDecay()
{
  // decay process
  G4Decay* decay = new G4Decay();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() )
  {
    G4ParticleDefinition* particle = particleIterator->value();
    if (decay->IsApplicable(*particle) && !particle->IsShortLived())
    {
      ph->RegisterProcess(decay, particle);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XraySPOPhysicsList::AddRadioactiveDecay()
{
  G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();

  radioactiveDecay->SetARM(true);                //Atomic Rearangement

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

