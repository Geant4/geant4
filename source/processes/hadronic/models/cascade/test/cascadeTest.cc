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
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4Proton.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4BertiniCascade.hh"
#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"

int main() {
 
  // Set test parameters
  G4int verboseLevel                = 2;                          
  G4int numberOfCascades            = 10; 
  G4double projectileMomentum       = 1.0 * GeV;
  //G4ParticleDefinition *particle = G4PionMinus::PionMinus();  
  //G4ParticleDefinition *particle = G4Neutron::Neutron();  
  G4ParticleDefinition *particle = G4Proton::Proton();  

  // Create projectile particle
  G4DynamicParticle *projectile = new G4DynamicParticle(); 
  projectile->SetDefinition(particle);
  //projectile->SetKineticEnergy( 1.0 * GeV);
  projectile->SetMomentum(projectileMomentum);
  projectile->SetMomentumDirection(1.0, 0.0, 0.0);  

  if (verboseLevel > 1) {
    G4cout << "projectile" << G4endl;
    G4cout << " type           : " << projectile->GetDefinition()        << G4endl;
    G4cout << " kinetic energy : " << projectile->GetKineticEnergy()     << G4endl;
    G4cout << " momentum       : " << projectile->GetMomentum()          << G4endl;
    G4cout << " p direction    : " << projectile->GetMomentumDirection() << G4endl;
  }

  // Set projectile particle track
  G4ThreeVector v;                                            
  v.setX(0.0 * fermi); 
  v.setY(0.0 * fermi); 
  v.setZ(0.0 * fermi);
  G4Track aTrack(projectile, 0, v);

  // Set target nucleus
  G4Nucleus targetNucleus;                                        
  G4double a(10);
  G4double z(10);
  targetNucleus.SetParameters(a, z);
  G4VParticleChange *cascadeParticles;

  if (verboseLevel > 1) {
    G4cout << "target" << G4endl;
    G4cout << " a              : " << a                              << G4endl;
    G4cout << " z              : " << z                              << G4endl;
    G4cout << " atomic mass    : " << targetNucleus.AtomicMass(a, z) << G4endl;
  }

  G4BertiniCascade *theCascade  = new G4BertiniCascade();
  for (G4int cascadeID =1 ; cascadeID <= numberOfCascades; cascadeID++) { 
    if (verboseLevel > 1) G4cout << "inc " << cascadeID << G4endl;
    cascadeParticles = theCascade->ApplyYourself(aTrack, targetNucleus);
  }

  delete projectile;
  delete theCascade;   
}
