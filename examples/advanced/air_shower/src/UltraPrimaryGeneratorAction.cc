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

  // electron
  //    G4ParticleDefinition* electron = particleTable->FindParticle(particleName="e-");
  // pi-minus
  //  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="pi-");
  // proton
  //  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="proton");
  // mu-minus
  //    G4ParticleDefinition*     muon = particleTable->FindParticle(particleName="mu-");

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
  posDist->SetRadius(10.0*cm);

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

}
  particleGun->GeneratePrimaryVertex(anEvent);
}


