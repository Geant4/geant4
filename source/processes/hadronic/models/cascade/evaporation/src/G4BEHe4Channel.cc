// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000

#include "G4BEHe4Channel.hh"
#include "G4Alpha.hh"

G4BEHe4Channel::G4BEHe4Channel()
{
  name = "He4";
    verboseLevel = 0;
    particleA = 4;
    particleZ = 2;
    rho = 0.70588;
    A=4;
    spin = 0;
}


G4BEHe4Channel::~G4BEHe4Channel()
{
}


G4DynamicParticle * G4BEHe4Channel::emit()
{
  G4double u , v , w;
  G4DynamicParticle * pParticle = new G4DynamicParticle;
  pParticle -> SetDefinition( G4Alpha::Alpha() );
  pParticle -> SetKineticEnergy( sampleKineticEnergy() ); 
  isotropicCosines( u , v , w );
  pParticle -> SetMomentumDirection( u , v , w );  
  return pParticle;
}


G4double G4BEHe4Channel::coulombFactor()
{
  return 0;
}


G4double G4BEHe4Channel::qmFactor()
{
  // Coefficient representing the quantum mechanical barrier
  // penetration, see Dostrovsky, Phys. Rev. 116, 1959.
  return qmFactorForAlpha();
}
