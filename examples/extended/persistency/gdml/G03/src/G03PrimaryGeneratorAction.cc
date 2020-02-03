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
/// \file persistency/gdml/G03/src/G03PrimaryGeneratorAction.cc
/// \brief Implementation of the G03PrimaryGeneratorAction class
//
//
//
// Class G03PrimaryGeneratorAction implementation
//
// ----------------------------------------------------------------------------

#include "G03PrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G03PrimaryGeneratorAction::G03PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), 
   fParticleGun(0), 
   fParticleTable(0)
{
  // Particle gun and particle table 
  //
  fParticleGun = new G4ParticleGun();
  fParticleTable = G4ParticleTable::GetParticleTable();

  // Default particle
  //
  fParticleGun->SetParticleDefinition(fParticleTable->FindParticle("geantino"));
  fParticleGun->SetParticleEnergy( 1.0*MeV );

  G4ThreeVector err1=G4ThreeVector(-1260,-560,40); // outside
  G4ThreeVector err2=G4ThreeVector(100,-240,120);  // inside
  G4ThreeVector err2v=(err2-err1).unit();
  
  fParticleGun->SetParticleMomentumDirection(err2v);
  fParticleGun->SetParticlePosition(err1);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G03PrimaryGeneratorAction::~G03PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
