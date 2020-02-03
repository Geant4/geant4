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
/// \file GB07/src/GB07PrimaryGeneratorAction.cc
/// \brief Implementation of the GB07PrimaryGeneratorAction class
//
<<<<<<< HEAD:examples/extended/field/field02/src/F02RunAction.cc
//
// $Id: F02RunAction.cc 92497 2015-09-02 07:23:12Z gcosmo $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F02RunAction.hh"
#include "F02RunMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
=======
#include "GB07PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/biasing/GB07/src/GB07PrimaryGeneratorAction.cc

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB07PrimaryGeneratorAction::GB07PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="proton");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(10.*GeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-50*cm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB07PrimaryGeneratorAction::~GB07PrimaryGeneratorAction()
{
<<<<<<< HEAD:examples/extended/field/field02/src/F02RunAction.cc
  // save Rndm status
  if (fSaveRndm > 0)
  {
      G4Random::showEngineStatus();
      G4Random::saveEngineStatus("beginOfRun.rndm");
  }
  G4UImanager* ui = G4UImanager::GetUIpointer();
 
  G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();

  if (visManager) ui->ApplyCommand("/vis/scene/notifyHandlers");
=======
  delete fParticleGun;
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/biasing/GB07/src/GB07PrimaryGeneratorAction.cc
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB07PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
<<<<<<< HEAD:examples/extended/field/field02/src/F02RunAction.cc
  if (G4VVisManager::GetConcreteInstance())
  {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  // save Rndm status

  if (fSaveRndm == 1)
  {
    G4Random::showEngineStatus();
    G4Random::saveEngineStatus("endOfRun.rndm");
  }
=======
  fParticleGun->GeneratePrimaryVertex(anEvent);
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/biasing/GB07/src/GB07PrimaryGeneratorAction.cc
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
