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
// $Id: Tst50PrimaryGeneratorAction.cc,v 1.4 2002-12-18 17:04:42 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst50PrimaryGeneratorAction.hh"
//#include "Tst50DetectorConstruction.hh"
#include "Tst50PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50PrimaryGeneratorAction::Tst50PrimaryGeneratorAction()

{
  G4int n_particle = 1;

  particleGun = new G4ParticleGun(n_particle);
  gunMessenger = new Tst50PrimaryGeneratorMessenger(this);
  
energy=3.*MeV;   
// default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
   G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50PrimaryGeneratorAction::~Tst50PrimaryGeneratorAction()
{
  delete gunMessenger;
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  
  particleGun->SetParticlePosition(G4ThreeVector(0.*m,0.*m,-1.*m));
  
  particleGun->GeneratePrimaryVertex(anEvent);
  

}









//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




