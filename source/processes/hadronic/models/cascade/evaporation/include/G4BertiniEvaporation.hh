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

#ifndef G4VEvaporation_h
#define G4VEvaporation_h 1

#include "globals.hh"
#include "G4LayeredNucleus.hh"
#include "G4BertiniEvaporationChannel.hh"
#include "G4ParticleChange.hh"   


class G4BertiniEvaporation 
{

public:
  G4BertiniEvaporation();
  ~G4BertiniEvaporation(); 

  G4VParticleChange * BreakItUp( G4LayeredNucleus & nucleus);
  void setVerboseLevel( const G4int verbose );
  
private:  
  G4int verboseLevel;
  vector< G4BertiniEvaporationChannel * > channelVector;
  void G4BertiniEvaporation::fillParticleChange( vector< G4DynamicParticle * > secondaryParticleVector,
						 G4ParticleChange * theParticleChange );
  void G4BertiniEvaporation::splitBe8( const G4double E, 
				       const G4ThreeVector boost,
				       vector< G4DynamicParticle * > & secondaryParticleVector);
  void G4BertiniEvaporation::isotropicCosines( G4double & u, G4double & v,G4double & w );
};

#endif
