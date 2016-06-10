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
// $Id: UltraActionInitializer.cc 66241 2012-12-13 18:34:42Z gunter $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "UltraActionInitializer.hh"

#include "UltraPrimaryGeneratorAction.hh"
#include "UltraRunAction.hh"
#include "UltraEventAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UltraActionInitializer::UltraActionInitializer() : 
  G4VUserActionInitialization()
{
  //Note: We instantiate a GPS here, because we want to set the defaults.
  //Since GPS distributions are shared among threads, we need to be sure
  //that defaults are done only once.
  //See comments in file G4GeneralParticleSource.hh
    masterGPS = new G4GeneralParticleSource();
    // Define here the user default properties for the General Particle Source (GPS)
     // Can be modified through the GPS Messenger (/gps/... commands)

     G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
     G4String particleName;

     G4ParticleDefinition* opticalphoton = particleTable->FindParticle(particleName="opticalphoton");

     //....PARTICLE DEFINITIONS
     masterGPS->SetParticleDefinition(opticalphoton);

     G4ThreeVector Polarization = G4ThreeVector(1.,1.,0.) ;
     masterGPS->SetParticlePolarization(Polarization);

     // DEFINE A MONO-ENERGETIC SOURCE
     G4SPSEneDistribution *eneDist = masterGPS->GetCurrentSource()->GetEneDist() ;
     eneDist->SetEnergyDisType("Mono");
     eneDist->SetMonoEnergy(3.0*eV);

     // SET POSITION DISTRIBUTION
     G4SPSPosDistribution *posDist = masterGPS->GetCurrentSource()->GetPosDist() ;
     posDist->SetPosDisType("Plane");
     posDist->SetPosDisShape("Circle");
     posDist->SetRadius(20.0*cm);

   #ifdef ULTRA_MIRROR_USE
   #define ULTRA_REFLECTION_USE
   #endif

   #ifdef ULTRA_GROUND_USE
   #define ULTRA_REFLECTION_USE
   #endif

     G4SPSAngDistribution *angDist = masterGPS->GetCurrentSource()->GetAngDist() ;

   #ifdef ULTRA_REFLECTION_USE
     angDist->SetParticleMomentumDirection(G4ThreeVector(0.0,-1.0,0.0)) ;
     posDist->SetPosRot1(G4ThreeVector(1.,0.,0.));
     posDist->SetPosRot2(G4ThreeVector(0.,0.,-1.));
     posDist->SetCentreCoords(G4ThreeVector(0.0*cm,90.0*cm,150.0*cm));

   #else
     angDist->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,-1.0)) ;
     posDist->SetCentreCoords(G4ThreeVector(0.0*cm,0.0*cm,150.0*cm));

   #endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UltraActionInitializer::~UltraActionInitializer() 
{
   delete masterGPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UltraActionInitializer::Build() const 
{
  
  // primary generator
  SetUserAction(new UltraPrimaryGeneratorAction());

  //Thread-local RunAction: same class, but code controlled by IsMaster()
  SetUserAction(new UltraRunAction());
  SetUserAction(new UltraEventAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UltraActionInitializer::BuildForMaster() const
{ 
  //Thread-local RunAction: same class, but code controlled by IsMaster()
  SetUserAction(new UltraRunAction());
}

