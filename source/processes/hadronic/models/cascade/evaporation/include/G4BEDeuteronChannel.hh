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

#ifndef G4BEDeuteronChannel_h
#define G4BEDeuteronChannel_h 1

#include "G4BEChargedChannel.hh"

class G4BEDeuteronChannel : public G4BEChargedChannel
{
public:
  G4BEDeuteronChannel();
  virtual ~G4BEDeuteronChannel(); 

  G4DynamicParticle * emit();
  G4double constant();
  virtual G4double coulombFactor();
  virtual G4double qmFactor( );
  
private:  
};


#endif
