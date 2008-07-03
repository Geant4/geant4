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
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraPrimaryGeneratorAction.cc
//    ****************************************************
//  
//    Class used in the definition of the optical photons source
//    A plane, circular source is used. Depending on the source position, optical
//    photons may reach the UVscope directly or after reflection. By default direct
//    incidence is used. The source parameters can be set directly in this class
//    or through the GeneralParticleSource  messenger class.
//
#include "UltraPrimaryGeneratorAction.hh"
#include "UltraDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraPrimaryGeneratorAction::UltraPrimaryGeneratorAction()
{

  particleGun = new G4GeneralParticleSource();


  // Define here the user default properties for the General Particle Source (GPS)
  // Can be modified through the GPS Messenger (/gps/... commands)

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  G4ParticleDefinition* opticalphoton = particleTable->FindParticle(particleName="opticalphoton");

  //....PARTICLE DEFINITIONS
  particleGun->SetParticleDefinition(opticalphoton);

  G4ThreeVector Polarization = G4ThreeVector(1.,1.,0.) ;
  particleGun->SetParticlePolarization(Polarization);
  
  // DEFINE A MONO-ENERGETIC SOURCE
  G4SPSEneDistribution *eneDist = particleGun->GetCurrentSource()->GetEneDist() ;
  eneDist->SetEnergyDisType("Mono");
  eneDist->SetMonoEnergy(3.0*eV);

  // SET POSITION DISTRIBUTION 
  G4SPSPosDistribution *posDist = particleGun->GetCurrentSource()->GetPosDist() ;
  posDist->SetPosDisType("Plane");
  posDist->SetPosDisShape("Circle");
  posDist->SetRadius(20.0*cm);

#ifdef ULTRA_MIRROR_USE
#define ULTRA_REFLECTION_USE
#endif

#ifdef ULTRA_GROUND_USE
#define ULTRA_REFLECTION_USE
#endif

  G4SPSAngDistribution *angDist = particleGun->GetCurrentSource()->GetAngDist() ;

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraPrimaryGeneratorAction::~UltraPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4int iEvent = anEvent->GetEventID() ;
  if ( iEvent == 0 ){

    G4cout << particleGun->GetParticleDefinition()->GetParticleName()           << " " ;
    G4cout << particleGun->GetCurrentSource()->GetEneDist()->GetEnergyDisType() << " " ;
    G4cout << particleGun->GetCurrentSource()->GetPosDist()->GetPosDisType()    << G4endl ;


// Check if optical photon wavelength is within limits set for material optical properties tables. 

   


}
  particleGun->GeneratePrimaryVertex(anEvent);

    if (particleGun->GetParticleDefinition()->GetParticleName() == "opticalphoton"){
     
     	const UltraDetectorConstruction * detector =  
     	dynamic_cast<const UltraDetectorConstruction *>((G4RunManager::GetRunManager())->GetUserDetectorConstruction()) ;

	G4double lambda_min = detector->GetLambdaMin() ;
	G4double lambda_max = detector->GetLambdaMax() ;

       	G4double energy = particleGun->GetParticleEnergy() ;

	if (h_Planck*c_light/energy > lambda_max || h_Planck*c_light/energy < lambda_min){
	       G4cerr << "Error ! Optical photon energy (" << energy/eV << " eV) out of limits set by material optical properties tables. \n" 
              << "Please check that photon wavelength is within the following interval: [" 
              << lambda_min/nm << "," 
              << lambda_max/nm << "] nm" 
              << ", i.e., ["
              << h_Planck*c_light/lambda_max/eV << ","
              << h_Planck*c_light/lambda_min/eV << "] eV"
              << G4endl ;
	
             G4Exception("") ;
	}
 }

}


