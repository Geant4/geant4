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
// Derived from 
//  https://twiki.cern.ch/twiki/bin/view/Geant4/QuickMigrationGuideForGeant4V10
// Courtesy of A. Dotti
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/src/MyPrimaryGeneratorActionFromFile.cc
/// \brief Implementation of the MyPrimaryGeneratorActionFromFile class

#include "MyPrimaryGeneratorActionFromFile.hh"
#include "MyFileReader.hh"
#include "DetectorConstruction.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4AutoLock.hh"
#include "G4ParticleGun.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"

namespace { G4Mutex myPrimGenMutex = G4MUTEX_INITIALIZER; }

MyFileReader* MyPrimaryGeneratorActionFromFile::fileReader = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyPrimaryGeneratorActionFromFile::MyPrimaryGeneratorActionFromFile()
:G4VUserPrimaryGeneratorAction(), G4VStateDependent(),
 fParticleGun(0),
 fDetector(0)
{
  fDetector =
      dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
          ->GetUserDetectorConstruction());
   
  G4AutoLock lock(&myPrimGenMutex);
  if( !fileReader ) fileReader = new MyFileReader(); 
  fParticleGun = new G4ParticleGun(1);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(10.*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyPrimaryGeneratorActionFromFile::~MyPrimaryGeneratorActionFromFile()
{
  G4AutoLock lock(&myPrimGenMutex);
  if( fileReader ) { delete fileReader; fileReader = 0; }

  G4StateManager::GetStateManager()->DeregisterDependent(this);
  if(fParticleGun) delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyPrimaryGeneratorActionFromFile::GeneratePrimaries(G4Event* anEvent)
{
  //Energy is read from file
  //
  G4double nrj = 0;
  if(fileReader)
  {
    G4AutoLock lock(&myPrimGenMutex);
    nrj = fileReader->GetAnEvent();
  }
  fParticleGun->SetParticleEnergy(nrj*eV);

  //
  G4double thickness = fDetector->GetCytoThickness();
  G4double radius = fDetector->GetNuclRadius();

  G4double rx=1*m;
  G4double ry=1*m;
  G4double rz=1*m;
  G4double myRadius = 0;

  do
  {
    rx = (2*G4UniformRand()-1)*(radius+thickness)*1.01;
    ry = (2*G4UniformRand()-1)*(radius+thickness)*1.01;
    rz = (2*G4UniformRand()-1)*(radius+thickness)*1.01;
    myRadius = std::sqrt(rx*rx+ry*ry+rz*rz);

  } while (myRadius>radius+thickness || myRadius<radius) ;

  fParticleGun->SetParticlePosition(G4ThreeVector(rx,ry,rz));
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  
  //G4cout << "---> EVENT ID=" << anEvent->GetEventID()+1 << G4endl;
  //G4cout << "---> energy/eV=" << nrj << G4endl;
  
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool MyPrimaryGeneratorActionFromFile::Notify(G4ApplicationState requestedState)
{
  if(requestedState == G4State_Idle)
  {
    if(fParticleGun != 0) return true;

    fParticleGun  = new G4ParticleGun(1);

    // Define default primary
    G4ParticleDefinition* particle
             = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(300*eV);
    fParticleGun->SetParticlePosition(G4ThreeVector());
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  }

  return true;
}
