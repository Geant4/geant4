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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyPrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "HadrontherapyAnalysisManager.hh"

HadrontherapyPrimaryGeneratorAction::HadrontherapyPrimaryGeneratorAction()
{
  // Define the messenger
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
  
  // Define primary particles: protons
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable -> FindParticle("proton");
  particleGun -> SetParticleDefinition(particle); 

  // Define the energy of primary particles:
  // gaussian distribution with mean energy = 62.0 *MeV
  // and sigma = 400.0 *keV
  G4double defaultMeanKineticEnergy = 62.0 *MeV;
  meanKineticEnergy = defaultMeanKineticEnergy;

  G4double defaultsigmaEnergy = 400.0 *keV;
  sigmaEnergy = defaultsigmaEnergy;
  
#ifdef G4ANALYSIS_USE_ROOT
  // Write these values into the analysis if needed. Have to be written separately on change.
  HadrontherapyAnalysisManager::GetInstance()->setBeamMetaData(meanKineticEnergy, sigmaEnergy);
#endif

  // Define the parameters of the initial position: 
  // the y, z coordinates have a gaussian distribution
  
  G4double defaultX0 = -3000.0 *mm;
  X0 = defaultX0;

  G4double defaultY0 = 0.0 *mm;  
  Y0 = defaultY0;

  G4double defaultZ0 = 0.0 *mm;  
  Z0 = defaultZ0;

  G4double defaultsigmaY = 1. *mm;  
  sigmaY = defaultsigmaY;

  G4double defaultsigmaZ = 1. *mm;  
  sigmaZ = defaultsigmaZ;

  // Define the parameters of the momentum of primary particles: 
  // The momentum along the y and z axis has a gaussian distribution
  G4double defaultsigmaMomentumY = 0.0;  
  sigmaMomentumY = defaultsigmaMomentumY;

  G4double defaultsigmaMomentumZ = 0.0;  
  sigmaMomentumZ = defaultsigmaMomentumZ;
}

void HadrontherapyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
#ifdef G4ANALYSIS_USE_ROOT
  // Increment the event counter
  HadrontherapyAnalysisManager::GetInstance()->startNewEvent();
#endif

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
  G4double momentumX = 1.0;
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
 
  particleGun -> SetParticleMomentumDirection( G4ThreeVector(momentumX,momentumY,momentumZ) );

  // Generate a primary particle
  particleGun -> GeneratePrimaryVertex( anEvent ); 
} 

void HadrontherapyPrimaryGeneratorAction::SetmeanKineticEnergy (G4double val )  
{
	meanKineticEnergy = val;
#ifdef G4ANALYSIS_USE_ROOT
  // Update the beam-data in the analysis manager
  HadrontherapyAnalysisManager::GetInstance()->setBeamMetaData(meanKineticEnergy, sigmaEnergy);
#endif

} 

void HadrontherapyPrimaryGeneratorAction::SetsigmaEnergy (G4double val )  
{ 
	sigmaEnergy = val;
#ifdef G4ANALYSIS_USE_ROOT
  // Update the sigmaenergy in the metadata.
  HadrontherapyAnalysisManager::GetInstance()->setBeamMetaData(meanKineticEnergy, sigmaEnergy);
#endif
}

void HadrontherapyPrimaryGeneratorAction::SetXposition (G4double val )  
{ X0 = val;}

void HadrontherapyPrimaryGeneratorAction::SetYposition (G4double val )  
{ Y0 = val;}

void HadrontherapyPrimaryGeneratorAction::SetZposition (G4double val )  
{ Z0 = val;}

void HadrontherapyPrimaryGeneratorAction::SetsigmaY (G4double val )  
{ sigmaY = val;}

void HadrontherapyPrimaryGeneratorAction::SetsigmaZ (G4double val )  
{ sigmaZ = val;}

void HadrontherapyPrimaryGeneratorAction::SetsigmaMomentumY (G4double val )  
{ sigmaMomentumY = val;}

void HadrontherapyPrimaryGeneratorAction::SetsigmaMomentumZ (G4double val )  
{ sigmaMomentumZ = val;}

G4double HadrontherapyPrimaryGeneratorAction::GetmeanKineticEnergy(void)
{ return meanKineticEnergy;}

