// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelPrimaryGeneratorAction.cc                *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of xray_telescope Physics list
// - Based on Chandra and XMM models by S Magni and F Lei
// 
//
// **********************************************************************

#include <stdio.h>

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "Randomize.hh"

#include "XrayTelPrimaryGeneratorAction.hh"
#include "XrayTelDetectorConstruction.hh"
#include "XrayTelPrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelPrimaryGeneratorAction::XrayTelPrimaryGeneratorAction(
							     XrayTelDetectorConstruction* XrayTelDC)
  :XrayTelDetector(XrayTelDC),rndmFlag("aperture")
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
  
  if ( rndmFlag == "aperture" ) {
    G4ThreeVector Position = GetRandomRingPosition (Rmin, Rmax);
    Position.rotateY(-90*degree);
    G4double tem = sin(Tmax);
    tem = tem*tem;
    G4ThreeVector Direction = GetRandomDirection (tem);
    Direction.rotateY(-90*degree);
    Position += G4ThreeVector ( 780.1*cm, 0.0*cm, 0.0*cm);
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XrayTelPrimaryGeneratorAction::GetRandomRingPosition( 
								   G4double minRad, G4double maxRad ) 
{
  G4ThreeVector position;
  G4double xx;
  G4double yy;
  G4double R = 0.;
  while ( R >maxRad || R < minRad){
    G4double rand1 = G4UniformRand();
    G4double rand2 = G4UniformRand();
    xx = -maxRad+2*maxRad*rand1;
    yy = -maxRad+2*maxRad*rand2;
    R = sqrt (xx*xx + yy*yy);
  }
  position.setX(xx);
  position.setY(yy);
  position.setZ(0.*m);
  return (position);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XrayTelPrimaryGeneratorAction::GetRandomDirection(G4double maxT)
{
  G4ThreeVector retval;

  G4double CosTheta;   
  G4double SinTheta;
  G4double Phi; 
  G4double SinPhi;
  G4double CosPhi;
  G4double rand;

  //	selecting a random direction
  rand = G4UniformRand();
  SinTheta = sqrt (maxT*rand);

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....














