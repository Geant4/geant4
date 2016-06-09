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
// $Id: PrimaryGeneratorAction.cc,v 1.3 2010-07-16 07:37:48 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PrimaryGeneratorAction.hh"
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
{
  // default particle kinematic
  //
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
    
  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
  particleGun->SetParticleDefinition(particle);
        
  particleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  
  action1 = new PrimaryGeneratorAction1(particleGun);
  action2 = new PrimaryGeneratorAction2(particleGun);
  action3 = new PrimaryGeneratorAction3(particleGun);
  action4 = new PrimaryGeneratorAction4(particleGun);
  
  selectedAction = 1;
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete action1;
  delete action2;
  delete action3;
  delete action4;
  delete particleGun;  
  delete gunMessenger;      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  switch(selectedAction)
  {
   case 1:
    action1->GeneratePrimaries(anEvent);
    break;
   case 2:
    action2->GeneratePrimaries(anEvent);
    break;
   case 3:
    action3->GeneratePrimaries(anEvent);
    break;
   case 4:
    action4->GeneratePrimaries(anEvent);
    break;    
   default:
    G4cerr << "Invalid generator action" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
