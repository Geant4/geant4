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
// $Id: AXPETPrimaryGeneratorAction.cc,v 1.1 2008-09-03 13:34:03 gcosmo Exp $
// ------------------------------------------------------------
// Geant4 class implementation file
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#include "AXPETPrimaryGeneratorAction.hh"
#include "AXPETPrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

AXPETPrimaryGeneratorAction::AXPETPrimaryGeneratorAction()
{

  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new AXPETPrimaryGeneratorMessenger(this);
  
  //default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e-");

  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleTime(0.0*ns);

  particleGun->SetParticlePosition(G4ThreeVector(0.0*cm,0.0*mm,15.0*cm));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  particleGun->SetParticleEnergy(20.0*keV);
}

AXPETPrimaryGeneratorAction::~AXPETPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

void AXPETPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4double sigma_xy = 10.*mm;//2.5;// ~6 mm (FWHM)
  G4double x = G4RandGauss::shoot(0.,sigma_xy);
  G4double y = G4RandGauss::shoot(0.,sigma_xy);
  G4double z = G4RandGauss::shoot(0.,sigma_xy);//15.0*mm;
  G4ThreeVector dir;
  dir=(G4ThreeVector(0,0,0.0*mm)-G4ThreeVector(x,y,z));
  dir=dir/dir.mag();
  particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
  particleGun->SetParticleMomentumDirection(dir);
  G4cout << x << '\t' << y << '\t' << z << " position [mm] of Primaries" << G4endl;
  G4cout << dir.x()<<'\t'<<dir.y()<<'\t'<<dir.z()<<'\t' <<" direction of Primaries" << G4endl;
  particleGun->GeneratePrimaryVertex(anEvent);
}

void AXPETPrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

void AXPETPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 if (particleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton")
   {
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
               "the particleGun is not an opticalphoton" << G4endl;
     return;
   }
     	       
 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = particleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton); 
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product; 
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 particleGun->SetParticlePolarization(polar);
}

