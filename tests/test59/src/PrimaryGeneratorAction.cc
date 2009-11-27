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
// $Id: PrimaryGeneratorAction.cc,v 1.1 2009-11-27 16:06:28 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC) : Detector(DC)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  SetDefaultKinematic();
  rndmBeam = 0.;
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;  
}
  
void PrimaryGeneratorAction::SetDefaultKinematic()
{    
  // default particle kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(30.*MeV);
  G4double x0 = -0.5*(Detector->GetWorldSizeX());
  particleGun->SetParticlePosition(G4ThreeVector(x0, 0*cm, 0*cm));  
}

// this function is called at the begining of event. It randomizes the beam.
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  static const G4double tlim = 6.5*mrad;
  static const G4double sigm = 2.4*mrad;
  G4double theta = 1.;
  while (std::fabs(theta) > tlim) theta = CLHEP::RandGauss::shoot(0, sigm);
  G4double phi = twopi*G4UniformRand();
  particleGun->SetParticleMomentumDirection(G4ThreeVector(std::cos(theta),
                                                          std::sin(theta)*std::cos(phi),
                                                          std::sin(theta)*std::sin(phi)));
  if ( rndmBeam > 0.) 
  {
      G4ThreeVector oldPosition = particleGun->GetParticlePosition();    
      G4double rbeam = 0.5*(Detector->GetAbsorberSizeYZ())*rndmBeam;
      G4double x0 = oldPosition.x();
      G4double y0 = oldPosition.y() + (2*G4UniformRand()-1.)*rbeam;
      G4double z0 = oldPosition.z() + (2*G4UniformRand()-1.)*rbeam;
      particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
      particleGun->GeneratePrimaryVertex(anEvent);
      particleGun->SetParticlePosition(oldPosition);      
  }
  else  particleGun->GeneratePrimaryVertex(anEvent); 
}
