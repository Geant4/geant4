#define hTestPrimaryGeneratorAction_CPP 

//---------------------------------------------------------------------------
//
// ClassName:   hTestPrimaryGeneratorAction
//  
// Description: Generate primary beam 
//
// Authors:    0.6.04.01 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "hTestPrimaryGeneratorAction.hh"
#include "hTestPrimaryGeneratorMessenger.hh"
#include "Randomize.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorAction::hTestPrimaryGeneratorAction(
			     hTestDetectorConstruction* det):
  theDet(det)
{
  InitializeMe();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::InitializeMe()
{
  verbose = theDet->GetVerbose();
  theMessenger = new hTestPrimaryGeneratorMessenger(this);
  particleGun = new G4ParticleGun();
  counter = 0;
  x0 = 0.0; 
  y0 = 0.0;
  z0 = 0.0;
  sigmaX = 0.0;
  sigmaY = 0.0;
  sigmaZ = 0.0;
  sigmaE = 0.0;
  minCosTheta = 1.0;
  energy = 10.0*MeV;
  position  = G4ThreeVector(x0,y0,z0);
  direction = G4ThreeVector(0.0,0.0,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorAction::~hTestPrimaryGeneratorAction()
{
  delete particleGun;
  delete theMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  counter++ ;
  verbose = theDet->GetVerbose();

  // Simulation of beam position
  G4double x = x0;
  G4double y = y0;
  G4double z = z0;
  if(0.0 < sigmaX) x += G4RandGauss::shoot(0.0,sigmaX);
  if(0.0 < sigmaY) y += G4RandGauss::shoot(0.0,sigmaY);
  if(0.0 < sigmaZ) z += G4RandGauss::shoot(0.0,sigmaZ);
  position  = G4ThreeVector(x,y,z);
  particleGun->SetParticlePosition(position);

  // Simulation of beam direction
  G4double ux = direction.x();
  G4double uy = direction.y();
  G4double uz = direction.z();

  // Beam particles are uniformly distributed over phi, cosTheta 
  if(1.0 > minCosTheta) {
    uz = minCosTheta + (1.0 - minCosTheta)*G4UniformRand() ;
    ux = sqrt(1.0 - uz*uz) ;
    uy = ux ;
    G4double phi = 360.0*deg*G4UniformRand() ;
    ux *= cos(phi) ;
    uy *= sin(phi) ;
    direction = G4ThreeVector(ux,uy,uz) ;
  }
 
  particleGun->SetParticleMomentumDirection(direction.unit());

  // Simulation of beam kinetic energy
  energy = particleGun->GetParticleEnergy();
  if(0.0 < sigmaE) {
    energy += G4RandGauss::shoot(0.0,sigmaE);
    if(0.0 > energy) energy = 0.0;
    particleGun->SetParticleEnergy(energy);
  }  

  G4ParticleDefinition* particle = particleGun->GetParticleDefinition();
  G4String particleName = particle->GetParticleName() ;

  if(verbose > 0) {
    G4cout << "Event#  " << counter 
           << "  Beam particle is generated by hTestPrimaryGeneratorAction " 
           << G4endl;
    G4cout << "ParticleName= " << particleName 
           << "  PDGcode= " << particle->GetPDGEncoding()
           << G4std::setprecision(5) 
	   << "   KinEnergy(GeV)= "
	   << energy/GeV 
	   << "   x(mm)= "
	   << x/mm 
	   << " y(mm)= "
	   << y/mm 
	   << " z(mm)= "
	   << z/mm 
           << "   ux= " 
	   << ux 
	   << " uy= "
	   << uy
	   << " uz= "
	   << uz 
	   << endl;
    }

  particleGun->GeneratePrimaryVertex(anEvent);
  if(verbose > 1) G4cout << "hTestPrimaryGeneratorAction: BeamOn" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











