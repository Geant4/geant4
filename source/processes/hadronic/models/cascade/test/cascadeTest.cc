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
 
  G4int verboseLevel                = 2;                          // set test parameters
  G4int numberOfCascades            = 10; 
  G4double projectileMomentum       = 1.0 * GeV;
  //G4ParticleDefinition *particle = G4PionMinus::PionMinus();  
  //G4ParticleDefinition *particle = G4Neutron::Neutron();  
  G4ParticleDefinition *particle = G4Proton::Proton();  


  G4DynamicParticle *projectile = new G4DynamicParticle(); // create projectile particle
  projectile->SetDefinition(particle);
  // projectile->SetKineticEnergy( 1.0 * GeV);
  projectile->SetMomentum(projectileMomentum);
  projectile->SetMomentumDirection(1.0, 0.0, 0.0);  

  if (verboseLevel > 1) {
    G4cout << "projectile" << G4endl;
    G4cout << " type           : " << projectile->GetDefinition()        << G4endl;
    G4cout << " kinetic energy : " << projectile->GetKineticEnergy()     << G4endl;
    G4cout << " momentum       : " << projectile->GetMomentum()          << G4endl;
    G4cout << " p direction    : " << projectile->GetMomentumDirection() << G4endl;
  }

  G4ThreeVector v;                                             // set projectile particle track
  v.setX(0.0 * fermi); 
  v.setY(0.0 * fermi); 
  v.setZ(0.0 * fermi);
  G4Track aTrack(projectile, 0, v);

  G4Nucleus targetNucleus;                                        // set target nucleus
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
    















