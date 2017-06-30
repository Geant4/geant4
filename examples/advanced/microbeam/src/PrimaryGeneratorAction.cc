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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fParticleGun  = new G4ParticleGun(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4int numEvent;
  numEvent=anEvent->GetEventID()+1;
  G4double x0,y0,z0,theta,phi,xMom0,yMom0,zMom0,e0;

  // INITIAL BEAM POSITION
  
  z0=-10000*mm;
  x0=10*mm;
  y0=10*mm;

  G4double sizeMax = 0.5*micrometer; // INITIAL BEAM POSITION UNIFORMLY SPREAD ON A DISK
  while (! (std::sqrt(x0*x0+y0*y0)<= sizeMax) )
  {
    x0 = CLHEP::RandFlat::shoot(-sizeMax,sizeMax);
    y0 = CLHEP::RandFlat::shoot(-sizeMax,sizeMax);
  }
  
  // INITIAL BEAM ENERGY
  
  e0= G4RandGauss::shoot(3*MeV,5.0955e-5*MeV); // AIFIRA ENERGY RESOLUTION

  // INITIAL BEAM DIVERGENCE

  do {
       theta=std::acos(1-G4UniformRand()*(1.0-std::cos(1.1e-6)))*rad;
     } 
  while(theta>1.1e-6*rad);

  phi=CLHEP::twopi*G4UniformRand()*rad;

  xMom0=std::sin(theta)*std::cos(phi);
  yMom0=std::sin(theta)*std::sin(phi);
  zMom0=std::cos(theta);

  // VERBOSE
  
  G4cout 
  << "-> Event # " << numEvent 
  << " generated " 
  << G4endl;
  
/*
  G4cout 
  << "-> Event # " << numEvent 
  << " : THETA from Z axis (mrad) = " << theta*1000 
  << " -- PHI (deg) = " << phi*180/CLHEP::pi 
  << " -- x0 (um) = " << x0/micrometer 
  << " -- y0 (um) = " << y0/micrometer 
  << " -- z0 (m) = " << z0/m 
  << " -- e0 (MeV) = " << e0/MeV 
  << G4endl;
*/

  fParticleGun->SetParticleEnergy(e0);

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xMom0,yMom0,zMom0));

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  
  G4ParticleDefinition* particle=
    G4ParticleTable::GetParticleTable()->FindParticle("alpha");
  
  fParticleGun->SetParticleDefinition(particle);
  
  fParticleGun->GeneratePrimaryVertex(anEvent);

}
