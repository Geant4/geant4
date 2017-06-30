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
/// \file exoticphysics/dmparticle/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// $Id$

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "TestParameters.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
  :G4VUserPrimaryGeneratorAction(),fDetector(det)
{
  fParticleGun  = new G4ParticleGun(1);
  fParticleGun->SetParticleEnergy(100 * GeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // this function is called at the begining of event
  if( 0 == anEvent->GetEventID() ) 
  {
    G4double z0 = -0.5*(fDetector->GetWorldSizeZ()) + 1*um;
    G4cout<<"Event = "<<anEvent->GetEventID()<<"; z0 = "<<z0;
    fParticleGun->SetParticlePosition( G4ThreeVector( 0.0, 0.0, z0) );
    G4cout<<"  Primary direction = "
          <<fParticleGun->GetParticleMomentumDirection()<<G4endl;
    TestParameters* param = TestParameters::GetPointer();
    param->SetBeamParticle(fParticleGun->GetParticleDefinition());
    param->SetBeamEnergy(fParticleGun->GetParticleEnergy());
    param->SetPositionZ(z0);
  }  
  fParticleGun->GeneratePrimaryVertex( anEvent );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

