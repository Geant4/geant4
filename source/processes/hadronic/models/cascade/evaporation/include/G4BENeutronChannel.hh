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

#ifndef G4BENeutronChannel_h
#define G4BENeutronChannel_h 1

#include "globals.hh"
#include "G4BertiniEvaporationChannel.hh"

class G4BENeutronChannel : public G4BertiniEvaporationChannel
{
public:
  G4BENeutronChannel();
  virtual ~G4BENeutronChannel(); 

  virtual void calculateProbability();

  G4DynamicParticle * emit();
  G4double sampleKineticEnergy();

private:  
  G4double alpha();
  G4double beta();
};

#endif
