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

#include "G4BEHe3Channel.hh"
#include "G4He3.hh"

G4BEHe3Channel::G4BEHe3Channel()
{
  name = "He3";
    verboseLevel = 0;
    particleA = 3;
    particleZ = 2;
    rho = 0.70588;
    A=3;
    spin = 0.5;
}


G4BEHe3Channel::~G4BEHe3Channel()
{
}


G4DynamicParticle*  G4BEHe3Channel::emit()
{
  G4double u , v , w;
  G4DynamicParticle * pParticle = new G4DynamicParticle;
  pParticle -> SetDefinition( G4He3::He3() );
  pParticle -> SetKineticEnergy( sampleKineticEnergy() ); 
  isotropicCosines( u, v, w);
  pParticle -> SetMomentumDirection( u, v, w);  
  return pParticle;
}


G4double G4BEHe3Channel::coulombFactor()
{
  return 0;
}


G4double G4BEHe3Channel::qmFactor()
{
  // Coefficient representing the quantum mechanical barrier
  // penetration, see Dostrovsky, Phys. Rev. 116, 1959.
  return qmFactorForAlpha() - 0.06;
}
