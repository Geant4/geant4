//  XrayTelPrimaryGeneratorAction.cc
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <stdio.h>

#include "XrayTelPrimaryGeneratorAction.hh"
#include "XrayTelDetectorConstruction.hh"
#include "XrayTelPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelPrimaryGeneratorAction::XrayTelPrimaryGeneratorAction(
                               XrayTelDetectorConstruction* XrayTelDC)
                              :XrayTelDetector(XrayTelDC),rndmFlag("isotropic")
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new XrayTelPrimaryGeneratorMessenger(this);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);

  // Default parameters for random generation

  Rmin = 30.5*cm;
  Rmax = 35.5*cm;
  Tmax = 1*degree;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelPrimaryGeneratorAction::~XrayTelPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void XrayTelPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4String ParticleName = particleGun->GetParticleDefinition()->GetParticleName();
  //  cout << "Particle selected : " << ParticleName << endl;
  
  if ( rndmFlag == "isotropic" ) {
    //    cout << "Isotropic generator selected with radius : " << Rmin << " -  " << Rmax << endl;
    G4ThreeVector Position = GetRandomShellPosition (Rmin, Rmax);
    G4ThreeVector Direction = GetRandomDirection ();
    G4ThreeVector Radius = G4ThreeVector ( G4ThreeVector ( 0,0,-6 ) - Position);
    if ( (Direction * Radius) < 0 ) Direction *= -1.0;  
    particleGun->SetParticleMomentumDirection( Direction );
    particleGun->SetParticlePosition( Position );
  }
  else if ( rndmFlag == "beam" ) {
    //    cout << "Beam generator selected with radius : " << Rmax  << endl;
    // get only z position, x and y are randomized on a circle
    G4ThreeVector position = particleGun->GetParticlePosition();
    G4ThreeVector rnd = GetRandomPositionOnaCircle(Rmax);
    position.setX (rnd.x());
    position.setY (rnd.y());
    particleGun->SetParticlePosition ( position );
    // set particle direction along -z axes with small randomization 
    G4double rndangx = G4UniformRand() * 0.008;
    G4double rndangy = G4UniformRand() * 0.008;
    particleGun->SetParticleMomentumDirection ( G4ThreeVector(rndangx, rndangy, 1) );
  } 
  else if ( rndmFlag == "aperture" ) {
    //    cout << "Aperture generator selected with radius : " << Rmin << " -  " << Rmax << endl;
    G4ThreeVector Position = GetRandomRingPosition (Rmin, Rmax);
    Position.rotateY(-90*degree);
    G4double tem = sin(Tmax);
    tem = tem*tem;
    G4ThreeVector Direction = GetRandomDirection (tem);
    Direction.rotateY(-90*degree);
    Position += G4ThreeVector ( 730.1*cm, 0.0*cm, 0.0*cm);
    particleGun->SetParticleMomentumDirection( Direction );
    particleGun->SetParticlePosition( Position );
  }
  else if (rndmFlag == "point") {
    // do nothing here use the defaults
  }
  else {
    G4Exception ( "No primary generator type selected" );
  }

  particleGun->GeneratePrimaryVertex(anEvent);
}
 
 G4double XrayTelPrimaryGeneratorAction::GetRandomEnergy ( G4int ParticleCode )
{
  // ParticleCode = 1 --> electron
  // ParticleCode = 2 --> proton
  G4double tmp;
  if ( ParticleCode == 1 ) {
    tmp=ElectronRandomEnergy->shoot()*(7-0.04)+0.04;
    return tmp;
  }
  else if ( ParticleCode == 2 ) {
    tmp=ProtonRandomEnergy->shoot()*(300-0.1)+0.1;
    return tmp;
  } 
}

 G4ThreeVector XrayTelPrimaryGeneratorAction::GetRandomDirection() 
{
  G4ThreeVector retval;

  G4double CosTheta;   
  G4double SinTheta;
  
  G4double Phi; 
  G4double SinPhi;
  G4double CosPhi;
  
  G4double rand;

  //
  //	selecting a random direction
  //

  rand = G4UniformRand();

  CosTheta = 2.0*rand -1.0;
  SinTheta = sqrt (1.-CosTheta*CosTheta);
  rand = G4UniformRand();
  Phi = twopi*rand;
  SinPhi = sin (Phi);
  CosPhi = cos (Phi);
  retval.setX(SinTheta*CosPhi);
  retval.setY(SinTheta*SinPhi);
  retval.setZ(CosTheta);

  return retval;
}

 G4ThreeVector XrayTelPrimaryGeneratorAction::GetRandomDirection(G4double Tmax)
{
  G4ThreeVector retval;

  G4double CosTheta;   
  G4double SinTheta;
  
  G4double Phi; 
  G4double SinPhi;
  G4double CosPhi;
  
  G4double rand;

  //
  //	selecting a random direction
  //
  //  SinTheta = 1.;
  //while (SinTheta >= Tmax ) {
    rand = G4UniformRand();
    SinTheta = sqrt (Tmax*rand);
    //}
  CosTheta = sqrt (1.-SinTheta*SinTheta);
  rand = G4UniformRand();
  Phi = twopi*rand;
  SinPhi = sin (Phi);
  CosPhi = cos (Phi);
  retval.setX(SinTheta*CosPhi);
  retval.setY(SinTheta*SinPhi);
  retval.setZ(CosTheta);

  return retval;
}


 G4ThreeVector XrayTelPrimaryGeneratorAction::GetRandomShellPosition(
                                              G4double Rmin, G4double Rmax )
{
  G4double rand = G4UniformRand();  
  G4double R = pow ( rand, 1.0/3.0 )*( Rmax - Rmin ) + Rmin;
  G4ThreeVector position = GetRandomDirection() *= R;
  position += G4ThreeVector ( 0, 0, -6*m );
  return (position);
}

 G4ThreeVector XrayTelPrimaryGeneratorAction::GetRandomRingPosition( 
                                              G4double Rmin, G4double Rmax ) 
{
  G4ThreeVector position;
  G4double xx;
  G4double yy;
  G4double R = 0.;
  while ( R >Rmax || R < Rmin){
    G4double rand1 = G4UniformRand();
    G4double rand2 = G4UniformRand();
    xx = -Rmax+2*Rmax*rand1;
    yy = -Rmax+2*Rmax*rand2;
    R = sqrt (xx*xx + yy*yy);
  }
  position.setX(xx);
  position.setY(yy);
  position.setZ(0.*m);
  return (position);
}


 G4ThreeVector XrayTelPrimaryGeneratorAction::Get2DRandomDirection() 
{
  G4ThreeVector retval;

  G4double CosTheta;   
  G4double SinTheta;
  
  G4double Phi; 
  G4double SinPhi;
  G4double CosPhi;
  
  G4double rand;

  //
  //	scelgo una direzione di propagazione casuale
  //

  rand = G4UniformRand();
  Phi = twopi*rand;
  SinPhi = sin (Phi);
  CosPhi = cos (Phi);
  retval.setX(CosPhi);
  retval.setY(SinPhi);
  retval.setZ(0);

  return retval;
}

 G4ThreeVector XrayTelPrimaryGeneratorAction::GetRandomPositionOnaCircle( G4double R )
{
  G4double rand = G4UniformRand();
  G4double RR = pow ( rand, 1.0/2.0 ) * R;
  
  return ( Get2DRandomDirection() *= RR );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


