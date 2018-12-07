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
/// \file runAndEvent/RE01/src/RE01PrimaryGeneratorAction.cc
/// \brief Implementation of the RE01PrimaryGeneratorAction class
//
//

#include "RE01PrimaryGeneratorAction.hh"
#include "RE01PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4HEPEvtInterface.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01PrimaryGeneratorAction::RE01PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    fHEPEvt(0), fParticleGun(0), fMessenger(0)
{
  const char* filename = "pythia_event.data";
  fHEPEvt = new G4HEPEvtInterface(filename);

  G4int n_particle = 1;
  G4ParticleGun* particleGun = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="mu+");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
  particleGun->SetParticleEnergy(100.*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  fParticleGun = particleGun;

  fMessenger = new RE01PrimaryGeneratorMessenger(this);
  fUseHEPEvt = true;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4AutoLock.hh"
namespace {
  G4Mutex RE01PrimGenDestrMutex = G4MUTEX_INITIALIZER;
  G4Mutex RE01PrimGenMutex = G4MUTEX_INITIALIZER;
}

RE01PrimaryGeneratorAction::~RE01PrimaryGeneratorAction()
{
  G4AutoLock lock(&RE01PrimGenDestrMutex);
  if(fHEPEvt) {delete fHEPEvt; fHEPEvt = 0;}
  delete fParticleGun;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(fUseHEPEvt)
  {
    G4AutoLock lock(&RE01PrimGenMutex);
    fHEPEvt->GeneratePrimaryVertex(anEvent);
    if(anEvent->GetNumberOfPrimaryVertex()==0)
    { // End-of-file
      G4Exception("RE01PrimaryGeneratorAction::GeneratePrimaries",
                  "RE01_0001",RunMustBeAborted,"End of file detected.");
    }
  }
  else
  { fParticleGun->GeneratePrimaryVertex(anEvent); }
}
