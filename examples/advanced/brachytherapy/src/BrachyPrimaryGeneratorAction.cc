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
//    ********************************************
//    *                                          *
//    *    BrachyPrimaryGeneratorAction.cc       *
//    *                                          *
//    ********************************************

#include "BrachyPrimaryGeneratorAction.hh"
#include "BrachyAnalysisManager.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include <math.h>
#include "G4RunManager.hh"

//....

BrachyPrimaryGeneratorAction::BrachyPrimaryGeneratorAction()
{
 
  G4int NumParticles = 1;

  m_pParticleGun = new G4ParticleGun(NumParticles);


 
  m_pParticleGun = new G4ParticleGun(NumParticles);


 

}

//....

BrachyPrimaryGeneratorAction::~BrachyPrimaryGeneratorAction()
{
  if(m_pParticleGun)
    delete m_pParticleGun;


}

//....

void BrachyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();

  G4ParticleTable* pParticleTable = G4ParticleTable::GetParticleTable();
  G4String ParticleName = "gamma";
  G4ParticleDefinition* pParticle = pParticleTable->FindParticle(ParticleName);
 

  m_pParticleGun->SetParticleDefinition(pParticle);

  //  Random generation of gamma source point inside the Iodium core
  G4double x,y,z;
  G4double radius= 0.30*mm;
 
 
  do{
    x = (G4UniformRand()-0.5)*(radius)/0.5;
    y = (G4UniformRand()-0.5)*(radius)/0.5;
  }while(x*x+y*y > radius*radius);
 
  z = (G4UniformRand()-0.5)*1.75*mm/0.5 ;

  G4ThreeVector position(x,y,z);
  m_pParticleGun->SetParticlePosition(position);


  // Random generation of the impulse direction
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
  m_pParticleGun->SetParticleMomentumDirection(direction);


  //Energy = 0.356*MeV;

  Energy = 356*keV;
  m_pParticleGun->SetParticleEnergy(Energy);

  //Check the energy
  analysis->Spectrum(Energy);
  
  m_pParticleGun->GeneratePrimaryVertex(anEvent);

 

}

