//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// $Id: HadrontherapyPositronPrimaryGeneratorAction.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyPrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

HadrontherapyPrimaryGeneratorAction::HadrontherapyPrimaryGeneratorAction()
{
  gunMessenger = new HadrontherapyPrimaryGeneratorMessenger(this);
  particleGun  = new G4ParticleGun();
  SetDefaultPrimaryParticle();  
}  

HadrontherapyPrimaryGeneratorAction::~HadrontherapyPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}
  
void HadrontherapyPrimaryGeneratorAction::SetDefaultPrimaryParticle()
{    
  // ****************************
  // Default primary particle
  // ****************************
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  particleGun->SetParticleDefinition(particle); 
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));

  G4double defaultMeanKineticEnergy = 63.450 *MeV;
  meanKineticEnergy = defaultMeanKineticEnergy;

  G4double defaultsigmaEnergy = 400.0 *keV;
  sigmaEnergy = defaultsigmaEnergy;

  G4double defaultX0 = -3248.59 *mm;  
  X0 = defaultX0;

  G4double defaultY0 = 0.0 *mm;  
  Y0 = defaultY0;

  G4double defaultZ0 = 0.0 *mm;  
  Z0 = defaultZ0;

  G4double defaultsigmaY = 1 *mm;  
  sigmaY = defaultsigmaY;

  G4double defaultsigmaZ = 1 *mm;  
  sigmaZ = defaultsigmaZ;

  G4double defaultsigmaMomentumY = 0.0001;  
  sigmaMomentumY = defaultsigmaMomentumY;

  G4double defaultsigmaMomentumZ = 0.0001;  
  sigmaMomentumZ = defaultsigmaMomentumZ;
}

void HadrontherapyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // ****************************************
  // Set the beam angular apread 
  // and spot size
  // beam spot size
  // ****************************************

  // Set the position of the primary particles
  G4double x = X0;
  G4double y = Y0;
  G4double z = Z0;

  if ( sigmaY > 0.0 )
    {
      y += G4RandGauss::shoot( Y0, sigmaY );
    }
 
  if ( sigmaZ > 0.0 )
    {
      z += G4RandGauss::shoot( Z0, sigmaZ );
    }

  particleGun -> SetParticlePosition(G4ThreeVector( x , y , z ) );
 
  // ********************************************
  // Set the beam energy and energy spread
  // ********************************************

  G4double kineticEnergy = G4RandGauss::shoot( meanKineticEnergy, sigmaEnergy );
  particleGun -> SetParticleEnergy ( kineticEnergy );

  // Set the direction of the primary particles
  G4double momentumY = 0.0;
  G4double momentumZ = 0.0;

  if ( sigmaMomentumY  > 0.0 )
    {
      momentumY += G4RandGauss::shoot( 0., sigmaMomentumY );
    }
  if ( sigmaMomentumZ  > 0.0 )
    {
      momentumZ += G4RandGauss::shoot( 0., sigmaMomentumZ );
    }
 
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1,momentumY,momentumZ));
  // Generate a primary particle
  particleGun -> GeneratePrimaryVertex( anEvent ); 
} 


void HadrontherapyPrimaryGeneratorAction::SetmeanKineticEnergy (G4double  val )  
{ meanKineticEnergy = val;} 

void HadrontherapyPrimaryGeneratorAction::SetsigmaEnergy (G4double  val )  
{ sigmaEnergy = val;}

void HadrontherapyPrimaryGeneratorAction::SetXposition (G4double  val )  
{ X0 = val;}

void HadrontherapyPrimaryGeneratorAction::SetYposition (G4double  val )  
{ Y0 = val;}

void HadrontherapyPrimaryGeneratorAction::SetZposition (G4double  val )  
{ Z0 = val;}

void HadrontherapyPrimaryGeneratorAction::SetsigmaY (G4double  val )  
{ sigmaY = val;}

void HadrontherapyPrimaryGeneratorAction::SetsigmaZ (G4double  val )  
{ sigmaZ = val;}

void HadrontherapyPrimaryGeneratorAction::SetsigmaMomentumY (G4double  val )  
{ sigmaMomentumY = val;}

void HadrontherapyPrimaryGeneratorAction::SetsigmaMomentumZ (G4double  val )  
{ sigmaMomentumZ = val;}
