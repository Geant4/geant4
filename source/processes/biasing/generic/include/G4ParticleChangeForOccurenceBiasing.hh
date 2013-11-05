#ifndef G4ParticleChangeForOccurenceBiasing_hh
#define G4ParticleChangeForOccurenceBiasing_hh 1

#include "G4VParticleChange.hh"

class G4ParticleChangeForOccurenceBiasing : public G4VParticleChange {
public:
  G4ParticleChangeForOccurenceBiasing(G4String name);
  ~G4ParticleChangeForOccurenceBiasing();

public:
  void     SetOccurenceWeightForNonInteraction(G4double w)       {fOccurenceWeightForNonInteraction = w;}
  G4double GetOccurenceWeightForNonInteraction()           const {return fOccurenceWeightForNonInteraction;}
  void     SetOccurenceWeightForInteraction(G4double w)          {fOccurenceWeightForInteraction = w;}
  G4double GetOccurenceWeightForInteraction()              const {return fOccurenceWeightForInteraction;}

public:
  // -- set a wrapped particle change AND USE IT TO UPDATE this occurence particle change state:
  void               SetWrappedParticleChange(G4VParticleChange* wpc);
  G4VParticleChange* GetWrappedParticleChange()                      const {return fWrappedParticleChange;}
public:
  // -- collect the secondaries from the wrapped particle change, apply weight correction, and clear wrapped particle change:
  void StealSecondaries();
  
public:
  // -- from base class G4VParticleChange:
  virtual G4Step* UpdateStepForAtRest   (G4Step* step);
  virtual G4Step* UpdateStepForAlongStep(G4Step* step);
  virtual G4Step* UpdateStepForPostStep (G4Step* step);

public:
  const G4String& GetName() const {return fName;}

private:
  G4String                                       fName;
  G4VParticleChange*            fWrappedParticleChange;
  G4double           fOccurenceWeightForNonInteraction;
  G4double              fOccurenceWeightForInteraction;

  
};

#endif
