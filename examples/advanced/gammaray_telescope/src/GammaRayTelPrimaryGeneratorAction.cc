// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelPrimaryGeneratorAction.cc,v 1.3 2000-11-24 16:57:00 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelPrimaryGeneratorAction  ------
//           by  G.Santin, F.Longo & R.Giannitrapani (13 nov 2000)
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelPrimaryGeneratorAction.hh"

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorAction::GammaRayTelPrimaryGeneratorAction
(GammaRayTelDetectorConstruction* GammaRayTelDC)
  :GammaRayTelDetector(GammaRayTelDC),rndmFlag("off"),
   nSourceType(0),nSpectrumType(0)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new GammaRayTelPrimaryGeneratorMessenger(this);
  
  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  particleGun->SetParticleEnergy(30.*MeV);
  G4double position = 0.5*(GammaRayTelDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorAction::~GammaRayTelPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  G4double z0 = 0.5*(GammaRayTelDetector->GetWorldSizeZ());
  G4double x0 = 0.*cm, y0 = 0.*cm;
  
  G4ThreeVector pos0;
  G4ThreeVector dir0;
  G4ThreeVector vertex0 = G4ThreeVector(x0,y0,z0);
  
  dir0 = G4ThreeVector(0.,0.,-1.);

  G4double theta, phi, y, f;
  G4double theta0,phi0;
  
  switch(nSourceType) {
  case 0:
    particleGun->SetParticlePosition(vertex0);
    particleGun->SetParticleMomentumDirection(dir0);
    break;
  case 1:
    // GS: Generate random position on the 4PIsphere to create a unif. distrib.
    // GS: on the sphere
    phi = G4UniformRand() * 2.0 * M_PI;
    do {
      y = G4UniformRand()*1.0;
      theta = G4UniformRand() * M_PI;
      f = sin(theta);
    } while (y > f);
    vertex0 = G4ThreeVector(1.,0.,0.);
    vertex0.setMag(dVertexRadius);
    vertex0.setTheta(theta);
    vertex0.setPhi(phi);
    particleGun->SetParticlePosition(vertex0);

    dir0 = G4ThreeVector(1.,0.,0.);
    do {
      phi = G4UniformRand() * 2.0 * M_PI;
      do {
	y = G4UniformRand()*1.0;
	theta = G4UniformRand() * M_PI;
	f = sin(theta);
      } while (y > f);
      dir0.setPhi(phi);
      dir0.setTheta(theta);
    } while (vertex0.dot(dir0) >= -0.7 * vertex0.mag());
    particleGun->SetParticleMomentumDirection((G4ParticleMomentum)dir0);

    break;
  case 2:
    // GS: Generate random position on the upper semi-sphere z>0 to create a unif. distrib.
    // GS: on a plane
    phi = G4UniformRand() * 2.0 * M_PI;
    do {
      y = G4UniformRand()*1.0;
      theta = G4UniformRand() * M_PI/2;
      f = sin(theta) * cos(theta);
    } while (y > f);
    vertex0 = G4ThreeVector(1.,0.,0.);
    
    G4double xy = GammaRayTelDetector->GetWorldSizeXY();
    G4double z = GammaRayTelDetector->GetWorldSizeZ();
    
    if (dVertexRadius > xy*0.5)
      { 
	G4cout << "vertexRadius too big " << G4endl;
	G4cout << "vertexRadius setted to " << xy*0.45 << G4endl;
	dVertexRadius = xy*0.45;
      }
    
    if (dVertexRadius > z*0.5)
      { 
	G4cout << "vertexRadius too high " << G4endl;
	G4cout << "vertexRadius setted to " << z*0.45 << G4endl;
	dVertexRadius = z*0.45;
      }
    

    vertex0.setMag(dVertexRadius);
    vertex0.setTheta(theta);
    vertex0.setPhi(phi);
    
    // GS: Get the user defined direction for the primaries and
    // GS: Rotate the random position according to the user defined direction for the particle

    dir0 = particleGun->GetParticleMomentumDirection();
    if (dir0.mag() > 0.001) 
      {
	theta0 = dir0.theta();
	phi0   = dir0.phi();   
      }
    
    if (theta0!=0.) 
      {
	G4ThreeVector rotationAxis(1.,0.,0.);
	rotationAxis.setPhi(phi0+M_PI/2.);
	vertex0.rotate(theta0+M_PI,rotationAxis);
      }
    particleGun->SetParticlePosition(vertex0);
    break;
  }

  
  G4double pEnergy;
  
  switch(nSpectrumType) {
  case 0:
    break;
  case 1:
    break;
  case 2:
    do {
      y = G4UniformRand()*100000.0;
      pEnergy = G4UniformRand() * 10. * GeV;
      f = pow(pEnergy * (1/GeV), -4.);
    } while (y > f);
    
    particleGun->SetParticleEnergy(pEnergy);
    
    break;
  case 3:
    break;
  }

  particleGun->GeneratePrimaryVertex(anEvent);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....








