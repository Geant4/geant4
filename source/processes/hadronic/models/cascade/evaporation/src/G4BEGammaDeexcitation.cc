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

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4Nucleus.hh"  
#include "G4BEGammaDeexcitation.hh"


G4BEGammaDeexcitation::G4BEGammaDeexcitation()
{
}


G4BEGammaDeexcitation::~G4BEGammaDeexcitation()
{
}


void G4BEGammaDeexcitation::setVerboseLevel( G4int level )
{
  verboseLevel = level ;
}


void G4BEGammaDeexcitation::setNucleusA( G4int a )
{
  nucleusA = a;
}


void G4BEGammaDeexcitation::setNucleusZ( G4int z )
{
  nucleusZ = z;
}


void G4BEGammaDeexcitation::setExcitationEnergy( G4double energy )
{
  excitationEnergy = energy;
}


G4double  G4BEGammaDeexcitation::sampleKineticEnergy()
{
  return excitationEnergy  * G4UniformRand(); 
}


G4DynamicParticle *  G4BEGammaDeexcitation::emit()
{
  // Isotropic distribution assumed to gammas
  G4double u, v, w;
  G4DynamicParticle * pParticle = new G4DynamicParticle;
  pParticle -> SetDefinition( G4Gamma::Gamma() );
  pParticle -> SetKineticEnergy( sampleKineticEnergy() ); 
  isotropicCosines( u, v, w );
  pParticle -> SetMomentumDirection( u, v, w );  
  return pParticle;
}


void G4BEGammaDeexcitation::isotropicCosines( G4double & u, 
					      G4double & v, 
					      G4double & w )
{
  // Samples isotropic random direction cosines.
  G4double CosTheta = 1.0 - 2.0 * G4UniformRand();
  G4double SinTheta = sqrt( 1.0 - CosTheta * CosTheta );
  G4double Phi = twopi * G4UniformRand();

  u = cos( Phi ) * SinTheta;
  v = cos( Phi ) * CosTheta,
  w = sin( Phi );

  return;
}
