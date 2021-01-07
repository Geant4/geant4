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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorAction0.hh"
#include "PrimaryGeneratorAction1.hh"
#include "PrimaryGeneratorAction2.hh"
#include "PrimaryGeneratorAction3.hh"
#include "PrimaryGeneratorAction4.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   fAction0(0),
   fAction1(0),
   fAction2(0),
   fAction3(0),
   fAction4(0),
   fSelectedAction(0),
   fGunMessenger(0)
{
  // default particle kinematic
  //
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
    
  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
  fParticleGun->SetParticleDefinition(particle);
        
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  
  fAction0 = new PrimaryGeneratorAction0(fParticleGun);
  fAction1 = new PrimaryGeneratorAction1(fParticleGun);
  fAction2 = new PrimaryGeneratorAction2(fParticleGun);
  fAction3 = new PrimaryGeneratorAction3(fParticleGun);
  fAction4 = new PrimaryGeneratorAction4(fParticleGun);
  
  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fAction0;
  delete fAction1;
  delete fAction2;
  delete fAction3;
  delete fAction4;
  delete fParticleGun;  
  delete fGunMessenger;      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  switch(fSelectedAction)
  {
   case 0:
    fAction0->GeneratePrimaries(anEvent);
    break;
   case 1:
    fAction1->GeneratePrimaries(anEvent);
    break;
   case 2:
    fAction2->GeneratePrimaries(anEvent);
    break;
   case 3:
    fAction3->GeneratePrimaries(anEvent);
    break;
   case 4:
    fAction4->GeneratePrimaries(anEvent);
    break;    
   default:
    G4cerr << "Invalid generator fAction" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
