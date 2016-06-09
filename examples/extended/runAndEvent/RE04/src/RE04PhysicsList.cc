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
/// \file runAndEvent/RE04/src/RE04PhysicsList.cc
/// \brief Implementation of the RE04PhysicsList class
//
// $Id: $
//
#include "RE04PhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "HadronPhysicsQGSP_BERT.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4SystemOfUnits.hh"    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04PhysicsList::RE04PhysicsList(G4String& parWorldName)
:fpWorldName(parWorldName)
{
  defaultCutValue = 0.01*mm;
  G4int ver = 1;
  SetVerboseLevel(ver);

  // EM Physics
  this->RegisterPhysics( new G4EmStandardPhysics(ver) );

  // Synchrotron Radiation & GN Physics
  this->RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
  this->RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
  this->RegisterPhysics( new G4HadronElasticPhysics(ver) );

  // Hadron Physics
  this->RegisterPhysics( new HadronPhysicsQGSP_BERT(ver));

  // Stopping Physics
  this->RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
  this->RegisterPhysics( new G4IonPhysics(ver));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04PhysicsList::~RE04PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04PhysicsList::ConstructProcess()
{
  AddTransportation();

  AddParallelWorldProcess();

  G4PhysConstVector::iterator itr;
  for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
    (*itr)->ConstructProcess();
  }
}

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ChargedGeantino.hh"

#include "G4ParallelWorldProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04PhysicsList::AddParallelWorldProcess()
{
  // Add parallel world process
  G4ParallelWorldProcess* theParallelWorldProcess
      = new G4ParallelWorldProcess("paraWorldProc");
  theParallelWorldProcess->SetParallelWorld(fpWorldName);
  theParallelWorldProcess->SetLayeredMaterialFlag();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    if(particle!=G4ChargedGeantino::Definition()) {
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddProcess(theParallelWorldProcess);
      if(theParallelWorldProcess->IsAtRestRequired(particle))
      {pmanager->SetProcessOrdering(theParallelWorldProcess, idxAtRest, 9999);}
      pmanager->SetProcessOrdering(theParallelWorldProcess, idxAlongStep, 1);
      pmanager->SetProcessOrdering(theParallelWorldProcess, idxPostStep, 9999);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04PhysicsList::SetCuts()
{
  SetCutsWithDefault();
}
