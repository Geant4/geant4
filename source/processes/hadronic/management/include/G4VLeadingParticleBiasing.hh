#ifndef G4VLeadingParticleBiasing_h
#define G4VLeadingParticleBiasing_h

class G4VParticleChange;

class G4VLeadingParticleBiasing 
{
  public:
  
  virtual G4VParticleChange * Bias(G4VParticleChange * result) = 0;
};

#endif
