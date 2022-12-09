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
//

#include "FlashPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"
#include "globals.hh"

FlashPrimaryGeneratorAction::FlashPrimaryGeneratorAction() {

  particleGun = new G4GeneralParticleSource();

  SetDefaultPrimaryParticle();
}

FlashPrimaryGeneratorAction::~FlashPrimaryGeneratorAction() {
  delete particleGun;
}

void FlashPrimaryGeneratorAction::SetDefaultPrimaryParticle() {
  
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *particle = particleTable->FindParticle("e-");
  particleGun->SetParticleDefinition(particle);


  G4double defaultX0 = -100 * cm; //the beam propagates in the x direction from x=-100cm

  X0 = defaultX0;

  G4double defaultY0 = 0.0 * mm;
  Y0 = defaultY0;

  G4double defaultZ0 = 0.0 * mm;
  Z0 = defaultZ0;

  G4double defaultsigmaY = 2.5*mm;

  sigmaY = defaultsigmaY;

  G4double defaultsigmaZ = 2.5*mm;

  sigmaZ = defaultsigmaZ;

  G4double defaultTheta = 1.5 * deg;

  Theta = defaultTheta;
}

void FlashPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

  //Energy is defined in 9MeVEF.mac and 7MeVEF.mac for the 9 MeV and 7 MeV spectrum respectively

  G4double x = X0;//the beam has a gaussian shape
  G4double y = Y0;
  G4double z = Z0;

  if (sigmaY > 0.0) {
    y += G4RandGauss::shoot(Y0, sigmaY);
  }

  if (sigmaZ > 0.0) {
    z += G4RandGauss::shoot(Z0, sigmaZ);
  }

  particleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(x, y, z));

      particleGun->GetCurrentSource()->GetAngDist()->SetAngDistType("beam2d"); //define angular divergence
        particleGun->GetCurrentSource()->GetAngDist()->DefineAngRefAxes("angref1", G4ThreeVector( 0, 0,1.));//define reference axis for divergence
        
        particleGun->GetCurrentSource()->GetAngDist()->SetBeamSigmaInAngX(Theta);//set planar divergence
        particleGun->GetCurrentSource()->GetAngDist()->SetBeamSigmaInAngY(Theta);//set planar divergence
        particleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(1,0,0));//beam propagates in the x direction

  particleGun->GeneratePrimaryVertex(anEvent);
}
