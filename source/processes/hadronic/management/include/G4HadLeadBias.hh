#ifndef G4HadLeadBias_h
#define G4HadLeadBias_h

#include "G4VLeadingParticleBiasing.hh"
#include "g4std/vector"
#include "G4VParticleChange.hh"

class G4HadLeadBias : public G4VLeadingParticleBiasing
{
  public:
  virtual G4VParticleChange * Bias(G4VParticleChange * result);
  virtual ~G4HadLeadBias() {};
};

#endif
