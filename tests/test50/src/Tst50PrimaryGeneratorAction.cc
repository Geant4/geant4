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
// $Id: Tst50PrimaryGeneratorAction.cc,v 1.8 2003-01-17 17:14:15 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4IonTable.hh"
#include "Tst50PrimaryGeneratorAction.hh"
//#include "Tst50DetectorConstruction.hh"
#include "Tst50PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50PrimaryGeneratorAction::Tst50PrimaryGeneratorAction():
  rndmDirection("off")
{
  G4int n_particle = 1;

  particleGun = new G4ParticleGun(n_particle);
  gunMessenger = new Tst50PrimaryGeneratorMessenger(this);
  
  
// default particle
 
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //G4ParticleDefinition* particle = particleTable->FindParticle("e-");
    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(1.*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*m,0.*m,-1.*m));
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
  
  if(rndmDirection=="on")
  {

G4double a,b,c;
 G4double n;
 do{
   a = (G4UniformRand()-0.5)/0.5;
   b = (G4UniformRand()-0.5)/0.5; 
   c = (G4UniformRand()-0.5)/0.5;
   n = a*a+b*b+c*c;
   }while(n > 1 || n == 0.0);
 n = sqrt(n);
 a /= n;
 b /= n;
 c /= n;

 G4ThreeVector direction(a,b,c);
 particleGun->SetParticleMomentumDirection(direction);



  }
 if(rndmDirection=="off")
  {


 particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));



  }

  particleGun->GeneratePrimaryVertex(anEvent);
  

}

G4double Tst50PrimaryGeneratorAction::GetInitialEnergy()
{
  G4double energy= particleGun->GetParticleEnergy(); 
  return energy;
}
/*
//added by albe for theta distribution 8/1/2003
G4ThreeVector Tst50PrimaryGeneratorAction::GetInitialMomentumDirection()
{
   G4ThreeVector MomentumDirection = particleGun->GetParticleMomentumDirection(); 
  return MomentumDirection;
}
*/

G4String Tst50PrimaryGeneratorAction::GetParticle()
{
  G4String name= particleGun->GetParticleDefinition()->GetParticleName();
  return name;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




