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
//
// $Id: G3toG4PrimaryGeneratorAction.cc,v 1.5 2006-06-29 17:20:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "Randomize.hh"
#include "G3toG4PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

G3toG4PrimaryGeneratorAction::G3toG4PrimaryGeneratorAction(){
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle 
    = particleTable->FindParticle(particleName="chargedgeantino");
  particleGun->SetParticleDefinition(particle);
  G4ThreeVector direction(0, 0, 1);
  particleGun->SetParticleMomentumDirection(direction.unit());
  particleGun->SetParticleEnergy(1.*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm, 0.*cm,0.*cm));
}

G3toG4PrimaryGeneratorAction::~G3toG4PrimaryGeneratorAction(){
  delete particleGun;
}

void 
G3toG4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  //G4ThreeVector direction = GetRandomDirection();
  //particleGun->SetParticleMomentumDirection( direction.unit() ) ;
  G4cout << ">>>>>>>> Primary direction: " 
         << particleGun->GetParticleMomentumDirection() << G4endl;
  particleGun->GeneratePrimaryVertex(anEvent);
}

G4ThreeVector 
G3toG4PrimaryGeneratorAction::GetRandomDirection() {

  G4ThreeVector retval;

  G4double CosTheta;
  G4double SinTheta;

  G4double Phi;
  G4double SinPhi;
  G4double CosPhi;

  G4double rand;

  rand = G4UniformRand();

  CosTheta = 2.0*rand -1.0;
  SinTheta = std::sqrt (1.-CosTheta*CosTheta);
  rand = G4UniformRand();
  Phi = twopi*rand;
  SinPhi = std::sin (Phi);
  CosPhi = std::cos (Phi);
  retval.setX(SinTheta*CosPhi);
  retval.setY(SinTheta*SinPhi);
  retval.setZ(CosTheta);

  return retval;
}

