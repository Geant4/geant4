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
/// \file runAndEvent/RE02/src/RE02PrimaryGeneratorAction.cc
/// \brief Implementation of the RE02PrimaryGeneratorAction class
//
//
//

#include "RE02PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
RE02PrimaryGeneratorAction::RE02PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

// default particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.));
  fParticleGun->SetParticleEnergy(150.0*MeV);
//
// default beam position
  G4double position = -200./2.*cm;
//
// Initial beam spot size in sigma.; This is not a part of ParticleGun.
  fSigmaPosition = 10.* mm;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm, 0.*cm, position));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
RE02PrimaryGeneratorAction::~RE02PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
void RE02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 

  G4ThreeVector position = fParticleGun->GetParticlePosition();
  G4double dx = (G4UniformRand()-0.5)*fSigmaPosition;
  G4double dy = (G4UniformRand()-0.5)*fSigmaPosition;
  position.setX(dx);
  position.setY(dy);
  fParticleGun->SetParticlePosition(position);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
