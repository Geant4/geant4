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

#ifndef G4VGammaDeExcCh_h
#define G4VGammaDeExcCh_h 1

#include "globals.hh"

class G4BEGammaDeexcitation
{
public:
  G4BEGammaDeexcitation();
  virtual ~G4BEGammaDeexcitation();

  void setVerboseLevel( G4int verbose ); 

  void setNucleusA( G4int inputA );
  void setNucleusZ( G4int inputZ );
  void setExcitationEnergy( G4double inputE );

  G4DynamicParticle * emit();

private:  
  G4double sampleKineticEnergy();
  G4int verboseLevel;
  G4int nucleusA;
  G4int nucleusZ;
  G4double excitationEnergy;
  void isotropicCosines( G4double&,
			 G4double&,
			 G4double& );
};


#endif
