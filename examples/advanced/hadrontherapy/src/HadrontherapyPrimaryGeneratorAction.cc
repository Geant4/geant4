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
//
// $Id: HadrontherapyPrimaryGeneratorAction.cc,v 1.0
//
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
// --------------------------------------------------------------

#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyDetectorConstruction.hh"
//#include "HadrontherapyPrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
// -------------------------------------------------------------------------
HadrontherapyPrimaryGeneratorAction::HadrontherapyPrimaryGeneratorAction(
                                            HadrontherapyDetectorConstruction* HadrontherapyDC)
  :HadrontherapyDetector(HadrontherapyDC)
  

{
  particleGun  = new G4ParticleGun();
  SetDefaultKinematic();
    
  //create a messenger for this class
  //gunMessenger = new HadrontherapyPrimaryGeneratorMessenger(this);
}

// -------------------------------------------------------------------------
HadrontherapyPrimaryGeneratorAction::~HadrontherapyPrimaryGeneratorAction()
{
  delete particleGun;
  //delete gunMessenger;  
}
  
// ------------------------------------------------------------------------
void HadrontherapyPrimaryGeneratorAction::SetDefaultKinematic()
{    
  // ****************************
  // default particle kinematic
  // ****************************
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  particleGun -> SetParticleDefinition(particle);
  particleGun -> SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun -> SetParticleEnergy(63.4 *MeV); 
  particleGun -> SetParticlePosition(G4ThreeVector( -200.0 *cm, 0.0 *cm,0.0 *cm ));
}

// -------------------------------------------------------------------------
void HadrontherapyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // ****************************************
  // Set the beam angular apread 
  // and spot size
  // beam spot size
  // ****************************************

  G4double X0 = -2000 *mm;
  G4double Y0 = 0.0 *cm;
  G4double Z0 = 0.0 *cm;  

  G4double sigmaY = 1.0 *mm ;
  G4double sigmaZ = 1.0 *mm ;

  G4double x = X0;
  G4double y = Y0;
  G4double z = Z0;

  if ( sigmaY > 0.0 )
    {
      y += G4RandGauss::shoot( 0., sigmaY );
    }
 
  if ( sigmaZ > 0.0 )
    {
      z += G4RandGauss::shoot( Z0, sigmaZ );
    }

  particleGun -> SetParticlePosition(G4ThreeVector( x , y , z ) );
 
  // ********************************************
  // Set the beam energy and energy spread
  // ********************************************

  G4double meanKineticEnergy = 63.4 *MeV;
  G4double sigmaK = 300.0 *keV;
  G4double kineticEnergy = G4RandGauss::shoot ( meanKineticEnergy, sigmaK );
 
  particleGun -> SetParticleEnergy ( kineticEnergy );
  particleGun -> GeneratePrimaryVertex( anEvent ); 

  G4double momentumY = 0.0;
  G4double momentumZ = 0.0;
  G4double sigmaMomentumY = 0.0001;
  G4double sigmaMomentumZ = 0.0001;

  if ( sigmaMomentumY  > 0.0 )
    {
      momentumY += G4RandGauss::shoot( 0., sigmaMomentumY );
    }
  if ( sigmaMomentumZ  > 0.0 )
    {
      momentumZ += G4RandGauss::shoot( 0., sigmaMomentumZ );
    }
 
  particleGun -> SetParticleMomentumDirection(G4ThreeVector(1.,momentumY,momentumZ));
}


