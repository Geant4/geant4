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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include "G4SystemOfUnits.hh"
#include "IORTPrimaryGeneratorAction.hh"
#include "IORTPrimaryGeneratorMessenger.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
  
IORTPrimaryGeneratorAction::IORTPrimaryGeneratorAction()
{
  // Define the messenger
  gunMessenger = new IORTPrimaryGeneratorMessenger(this);

  particleGun  = new G4ParticleGun();

  SetDefaultPrimaryParticle();  
}  

IORTPrimaryGeneratorAction::~IORTPrimaryGeneratorAction()
{
  delete particleGun;

  delete gunMessenger;
}
  
void IORTPrimaryGeneratorAction::SetDefaultPrimaryParticle()
{    
  // ****************************
  // Default primary particle
  // ****************************
  
  // Define primary particles: electrons // protons
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable -> FindParticle("e-");   // ("proton")
  particleGun -> SetParticleDefinition(particle); 

  // Define the energy of primary particles:
  // gaussian distribution with mean energy = 10.0 *MeV   
  // and sigma = 400.0 *keV
  G4double defaultMeanKineticEnergy = 10.0 *CLHEP::MeV;   
  meanKineticEnergy = defaultMeanKineticEnergy;

  G4double defaultsigmaEnergy = 100.0 *CLHEP::keV;   
  sigmaEnergy = defaultsigmaEnergy;
 
  // Define the parameters of the initial position: 
  // the y, z coordinates have a gaussian distribution 
  
  G4double defaultX0 = -862.817 *CLHEP::mm;                 
  X0 = defaultX0;

  G4double defaultY0 = 0.0 *CLHEP::mm;  
  Y0 = defaultY0;

  G4double defaultZ0 = 0.0 *CLHEP::mm;  
  Z0 = defaultZ0;

  G4double defaultsigmaY = 1. *CLHEP::mm;  
  sigmaY = defaultsigmaY;

  G4double defaultsigmaZ = 1. *CLHEP::mm;  
  sigmaZ = defaultsigmaZ;

  // Define the parameters of the momentum of primary particles: 
  // The momentum along the y and z axis has a gaussian distribution
 
/* 
  G4double defaultsigmaMomentumY = 0.0;  
  sigmaMomentumY = defaultsigmaMomentumY;

  G4double defaultsigmaMomentumZ = 0.0;  
  sigmaMomentumZ = defaultsigmaMomentumZ;
*/

  G4double defaultTheta = 6.0 *CLHEP::deg;  
  Theta = defaultTheta;

}

void IORTPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
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
  

  /*                          
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

  */
  
               
  G4double Mx;    
  G4double My;
  G4double Mz;
  G4double condizione;
  
while (true)  {

  //Mx =  CLHEP::RandFlat::shoot(0.9,1);
  //My =  CLHEP::RandFlat::shoot(-0.1,0.1);
  //Mz =  CLHEP::RandFlat::shoot(-0.1,0.1);

  Mx =  CLHEP::RandFlat::shoot(0.7,1);
  My =  CLHEP::RandFlat::shoot(-0.3,0.3); // ranges good for 0<Theta<20
  Mz =  CLHEP::RandFlat::shoot(-0.3,0.3);

  condizione = std::sqrt(Mx*Mx + My*My + Mz*Mz);

 
  if (condizione < 1)  {
    Mx = Mx/condizione;
    My = My/condizione;
    Mz = Mz/condizione;


    if (Mx > std::cos(Theta)) { 
      break;
        }
    }
}
  
 
  particleGun -> SetParticleMomentumDirection( G4ThreeVector(Mx,My,Mz) );
  

  // Generate a primary particle
  particleGun -> GeneratePrimaryVertex( anEvent ); 
} 

void IORTPrimaryGeneratorAction::SetmeanKineticEnergy (G4double val )  
{
  meanKineticEnergy = val;
  G4cout << "The mean Kinetic energy of the incident beam has been changed to (MeV):" 
         << meanKineticEnergy/MeV << G4endl; 
} 

void IORTPrimaryGeneratorAction::SetsigmaEnergy (G4double val )  
{ 
 sigmaEnergy = val;
 G4cout << "The sigma of the kinetic energy of the incident beam has been changed to (MeV):"
        << sigmaEnergy/MeV << G4endl; 
}

void IORTPrimaryGeneratorAction::SetXposition (G4double val )  
{ X0 = val;}

void IORTPrimaryGeneratorAction::SetYposition (G4double val )  
{ Y0 = val;}

void IORTPrimaryGeneratorAction::SetZposition (G4double val )  
{ Z0 = val;}

void IORTPrimaryGeneratorAction::SetsigmaY (G4double val )  
{ sigmaY = val;}

void IORTPrimaryGeneratorAction::SetsigmaZ (G4double val )  
{ sigmaZ = val;}

/*
void IORTPrimaryGeneratorAction::SetsigmaMomentumY (G4double val )  
{ sigmaMomentumY = val;}

void IORTPrimaryGeneratorAction::SetsigmaMomentumZ (G4double val )  
{ sigmaMomentumZ = val;}
*/

void IORTPrimaryGeneratorAction::SetTheta (G4double val ) 
{ Theta = val;}


G4double IORTPrimaryGeneratorAction::GetmeanKineticEnergy(void)
{ return meanKineticEnergy;}

