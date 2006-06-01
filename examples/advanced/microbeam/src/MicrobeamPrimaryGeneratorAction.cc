// -------------------------------------------------------------------
// $Id: MicrobeamPrimaryGeneratorAction.cc,v 1.3 2006-06-01 22:25:20 sincerti Exp $
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

#include "MicrobeamPrimaryGeneratorAction.hh"
#include "MicrobeamDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamPrimaryGeneratorAction::MicrobeamPrimaryGeneratorAction(MicrobeamDetectorConstruction* DC)
  :Detector(DC)
{
  particleGun  = new G4ParticleGun(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamPrimaryGeneratorAction::~MicrobeamPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
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
    x0 = RandFlat::shoot(-sizeMax,sizeMax);
    y0 = RandFlat::shoot(-sizeMax,sizeMax);
  }
  
  // INITIAL BEAM ENERGY
  
  e0= G4RandGauss::shoot(3*MeV,5.0955e-5*MeV); // AIFIRA ENERGY RESOLUTION

  // INITIAL BEAM DIVERGENCE

  do {
       theta=std::acos(1-G4UniformRand()*(1.0-std::cos(1.1e-6)))*rad;
     } 
  while(theta>1.1e-6*rad);

  phi=2*M_PI*G4UniformRand()*rad;

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
  << " -- PHI (deg) = " << phi*180/M_PI 
  << " -- x0 (um) = " << x0/micrometer 
  << " -- y0 (um) = " << y0/micrometer 
  << " -- z0 (m) = " << z0/m 
  << " -- e0 (MeV) = " << e0/MeV 
  << G4endl;
*/

  particleGun->SetParticleEnergy(e0);

  particleGun->SetParticleMomentumDirection(G4ThreeVector(xMom0,yMom0,zMom0));

  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  
  G4ParticleDefinition* particle=
    G4ParticleTable::GetParticleTable()->FindParticle("alpha");
  
  particleGun->SetParticleDefinition(particle);
  
  particleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


