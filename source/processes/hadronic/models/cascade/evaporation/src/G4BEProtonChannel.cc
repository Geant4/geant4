// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, 30-May-2000

#include "G4BEProtonChannel.hh"
#include "G4Proton.hh"

G4BEProtonChannel::G4BEProtonChannel()
{
  name = "proton";
    verboseLevel = 0;
    particleA = 1;
    particleZ = 1;
    rho = 0; 
    A=1;
    spin = 0.5;
}


G4BEProtonChannel::~G4BEProtonChannel()
{
}


G4DynamicParticle*  G4BEProtonChannel::emit()
{
  G4double u , v , w;
  G4DynamicParticle * pParticle = new G4DynamicParticle;
  pParticle -> SetDefinition( G4Proton::Proton() );
  pParticle -> SetKineticEnergy( sampleKineticEnergy() ); 
  isotropicCosines( u , v , w );
  pParticle -> SetMomentumDirection( u , v , w );  
  return pParticle;
}


G4double G4BEProtonChannel::coulombFactor()
{
  return coulombFactorForProton();
}


G4double G4BEProtonChannel::qmFactor()
{
  // Coefficient representing the quantum mechanical barrier
  // penetration, see Dostrovsky, Phys. Rev. 116, 1959.
  return qmFactorForProton();
}
