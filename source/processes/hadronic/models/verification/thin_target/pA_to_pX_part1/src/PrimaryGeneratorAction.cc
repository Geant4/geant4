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
// $Id: PrimaryGeneratorAction.cc,v 1.1 2003-05-27 13:44:49 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "PrimaryGeneratorAction.hh"

#include "TargetConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4DecayTable.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction(TargetConstruction* TC)
  :Target(TC),rndmFlag("off")
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

  // default particle kinematics

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="proton");

  /*
  ////////////
  //
  // Test code
  // 
  ////////////

  G4ParticleDefinition* test
                    = particleTable->FindParticle(particleName="pi+");
  G4cout << "Particle name: " << test->GetParticleName() << G4endl;
  G4DecayTable* theDecayTable = test->GetDecayTable();
  G4int nChannels = theDecayTable->entries();
  G4cout << " # of decay channels = " << nChannels << G4endl;
  G4int n;
  G4int d;
  for (n = nChannels -1; n >= 0; n--) {
    G4VDecayChannel* theChannel = theDecayTable->GetDecayChannel(n);
    G4int nDaughters = theChannel->GetNumberOfDaughters();
    G4cout << " nDaughters = " << nDaughters << G4endl;
    
    G4ParticleDefinition* aDaughter;
    for (d = 0; d < nDaughters; d++) {
      aDaughter = theChannel->GetDaughter(d);
      G4cout << " Pointer = " << aDaughter << G4endl;
      G4cout << "Name = " << theChannel->GetDaughterName(d) << G4endl;
      G4cout << " Daughter name: " << aDaughter->GetParticleName() << G4endl;
    }
  } 

  //
  // end of test code 
  //
  */

  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(50.*MeV);
  G4double position = -0.5*(Target->GetWorldLength());
  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,position));

}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the beginning of the event
  // 
  G4double x0 = 0.;
  G4double y0 = 0.;
  G4double z0 = -0.5*(Target->GetWorldLength());
  //  if (rndmFlag == "on") {
  //    y0 = (G4UniformRand()-0.5)*mm;
  //    z0 = (G4UniformRand()-0.5)*mm;
  //  } 
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  particleGun->GeneratePrimaryVertex(anEvent);
}
