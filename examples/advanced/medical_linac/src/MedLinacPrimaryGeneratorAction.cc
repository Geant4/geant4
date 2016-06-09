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
// $Id: MedLinacPrimaryGeneratorAction.cc,v 1.6 2005/12/07 14:18:41 guatelli Exp $
//
//
// Code developed by: M. Piergentili
//
#include "MedLinacPrimaryGeneratorAction.hh"
#include "MedLinacPrimaryGeneratorMessenger.hh"

#ifdef G4ANALYSIS_USE
#include "MedLinacAnalysisManager.hh"
#endif
#include "MedLinacDetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4IonTable.hh"

MedLinacPrimaryGeneratorAction::MedLinacPrimaryGeneratorAction()
  :pEnergy(99.0*MeV),sigma(12.7*MeV)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

//create a messenger for this class
  gunMessenger = new MedLinacPrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="e-"));
  particleGun->SetParticleEnergy(pEnergy);

  //**********************************************************************
  //**********************************************************************
  particleGun->SetParticlePosition(G4ThreeVector(0.0*cm, 0.0*cm, 123.0*cm));
  //particleGun->SetParticlePosition(G4ThreeVector(0.0*cm, 0.0*cm, 160.09855*cm));
  //**********************************************************************
  //**********************************************************************

}

MedLinacPrimaryGeneratorAction::~MedLinacPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

void MedLinacPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double energy;

  //#ifdef G4ANALYSIS_USE
  //MedLinacAnalysisManager* analysis = MedLinacAnalysisManager::getInstance();
  //#endif
 
  //BeamLine
   //G4ThreeVector v(0.0,0.0,-1.0);
  //particleGun->SetParticleMomentumDirection(v);


 //uniform distribution in solid angle (-1<cosTheta<-0.9999; 0<phi<2pi)
  //

  //G4double cosTheta = -0.0001*G4UniformRand()-0.9999;
  G4double cosTheta = CLHEP::RandGauss::shoot(-1.,0.00003);
  G4double phi = twopi * G4UniformRand();

    G4double sinTheta = sqrt(1. - cosTheta*cosTheta);
    G4double ux = sinTheta*cos(phi);
    G4double uy = sinTheta*sin(phi);
    G4double uz = cosTheta;

   particleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));

  energy=pEnergy;
  energy = CLHEP::RandGauss::shoot(pEnergy,sigma);
  particleGun->SetParticleEnergy(energy);

  //G4cout <<"the energy of the primaries is  "  << energy <<" "<<"MeV"<< G4endl;
  //G4cout <<"the sigma of the energy is  "  << sigma <<" "<<"MeV"<< G4endl;
  
  particleGun->GeneratePrimaryVertex(anEvent);
}



