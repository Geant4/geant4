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
/// \file VG01PrimaryGeneratorAction.cc
/// \brief Implementation of the VG01PrimaryGeneratorAction class
//
//  Authors: J. Apostolakis & S. Wenzel (CERN)  2018-2021
//
//  Started from FullCMS code by Mihaly Novak (CERN) 2017  

#include "VG01PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "G4RandomDirection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VG01PrimaryGeneratorAction::VG01PrimaryGeneratorAction( EGeneratorMode mode )
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   fMode(mode)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fParticleGun->SetParticleDefinition(
               particleTable->FindParticle(particleName="geantino"));
  fParticleGun->SetParticleEnergy(1.0*GeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0) ); // -1.2*m, 0.1, 0.1));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VG01PrimaryGeneratorAction::~VG01PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VG01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ThreeVector dir(1.0,0.0,0.0);

  // G4ThreeVector FixedMode[3] = 
  //  { G4ThreeVector(1,0,0), G4ThreeVector(1,0.1,0), G4ThreeVector( 1, 0, 0.1)};
  G4ThreeVector AxisModeVectors[3] = 
      { G4ThreeVector(1,0,0), G4ThreeVector(0,1,0), G4ThreeVector( 0, 0, 1.0)};

  if( fMode == kFreeMode ){
    // Do not fix the direction 
    //   -- UI commands to particle mode can change it!
    fParticleGun->GeneratePrimaryVertex(anEvent);
    return; 
    //----
  } 
  else if( fMode == kFixedMode )
  {
    G4int i = anEvent->GetEventID() % 3;
    dir= G4ThreeVector(1.0,0.0,0.0);
    switch(i)
    {
       case 0:
          break;
       case 1:
          dir.setY(0.1);
          break;
       case 2:
          dir.setZ(0.1);
          break;
    }
  } 
  else if( fMode == kUniformMode )
  {
    dir = G4RandomDirection();
  }
  else if ( fMode == kAxisMode )
  {
    G4int i = anEvent->GetEventID() % 3;
    dir = AxisModeVectors[i];
  }

  fParticleGun->SetParticleMomentumDirection(dir);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
