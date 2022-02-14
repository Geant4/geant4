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
/// \file medical/dna/range/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// $Id: PrimaryGeneratorAction.cc 73024 2013-08-15 09:11:40Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "G4RandomDirection.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
:G4VUserPrimaryGeneratorAction(), G4VStateDependent(),
 fParticleGun(0),
 fDetector(0)
{
  fDetector =
      dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
          ->GetUserDetectorConstruction());
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  G4StateManager::GetStateManager()->DeregisterDependent(this);
  if(fParticleGun) delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{    
  G4double R = fDetector->GetNPRadius();
  G4double Ryz=9999;

  //G4ParticleDefinition *particle = G4Electron::ElectronDefinition();
  //fParticleGun->SetParticleDefinition(particle);

  G4double scale =1.;

  G4double xpos=0,ypos=0,zpos=0;
  while (!(Ryz<R)){
    ypos = (2.*G4UniformRand()-1.) * R;
    zpos = (2.*G4UniformRand()-1.) * R;
    Ryz = std::sqrt(ypos*ypos+zpos*zpos);
  }
  xpos = std::sqrt(R*R-(ypos*ypos+zpos*zpos));
  fParticleGun->SetParticlePosition(G4ThreeVector(
                               -xpos*scale,ypos*scale,zpos*scale));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PrimaryGeneratorAction::Notify(G4ApplicationState requestedState)
{
  if(requestedState == G4State_Idle)
  {
    if(fParticleGun != 0) return true;
      G4double R = fDetector->GetNPRadius();
      G4double Ryz=9999;;

      //G4ParticleDefinition *particle = G4Electron::ElectronDefinition();
      //fParticleGun->SetParticleDefinition(particle);

      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-1,0,0));

      while (!(Ryz<R)){
        G4double ypos = (2.*G4UniformRand()-1.) * R;
        G4double zpos = (2.*G4UniformRand()-1.) * R;
        Ryz = std::sqrt(ypos*ypos+zpos*zpos);
        G4double xpos = std::sqrt(R*R-(ypos*ypos+zpos*zpos));
        fParticleGun->SetParticlePosition(G4ThreeVector(-xpos,ypos,zpos));
      }
      //fParticleGun->SetParticlePosition(G4ThreeVector(-Xposition,0,0));
      //fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  return true;
}

