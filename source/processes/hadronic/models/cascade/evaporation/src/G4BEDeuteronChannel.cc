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

#include "G4BEDeuteronChannel.hh"
#include "G4Deuteron.hh"


G4BEDeuteronChannel::G4BEDeuteronChannel()
{
  name = "deuteron";
  verboseLevel = 0;
  particleA = 2;
  particleZ = 1;
  rho = 0.70588;
  exmass = G4NucleiProperties::GetMassExcess( particleA, particleZ ); 
  A=2;
  spin = 1;
}


G4BEDeuteronChannel::~G4BEDeuteronChannel()
{
}


G4DynamicParticle*  G4BEDeuteronChannel::emit()
{
  G4double u, v, w;
  G4DynamicParticle * pParticle = new G4DynamicParticle;
  pParticle -> SetDefinition( G4Deuteron::Deuteron() );
  pParticle -> SetKineticEnergy( sampleKineticEnergy() ); 
  isotropicCosines( u, v, w );
  pParticle -> SetMomentumDirection( u, v, w );  
  return pParticle;
}


G4double G4BEDeuteronChannel::coulombFactor()
{
  // Coefficient c_i representing the variation of charged-particle
  // capture cross sections with Z for each particle.  See Dostrovsky,
  // Phys. Rev. 116, 1959.
  return 0.5 * coulombFactorForProton();
}


G4double G4BEDeuteronChannel::qmFactor()
{
  // Coefficient k_i representing the quantum mechanical barrier
  // penetration, see Dostrovsky, Phys. Rev. 116, 1959.
  return qmFactorForProton() + 0.06;
}
