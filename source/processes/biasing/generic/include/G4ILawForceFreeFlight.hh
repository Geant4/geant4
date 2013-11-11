#ifndef G4ILawForceFreeFlight_hh
#define G4ILawForceFreeFlight_hh 1

#include "G4VBiasingInteractionLaw.hh"

class G4ILawForceFreeFlight : public G4VBiasingInteractionLaw
{
public:
  G4ILawForceFreeFlight(G4String name = "forceFreeFlightLaw");
  virtual ~G4ILawForceFreeFlight();
  
public:
  virtual G4double     ComputeEffectiveCrossSectionAt(G4double               length) const;
  virtual G4double ComputeNonInteractionProbabilityAt(G4double               length) const;
  // -- sample the distribution
  virtual  G4double           SampleInteractionLength();
  // -- move by true path length, this position becomes the new initial point
  virtual G4double     UpdateInteractionLengthForStep(G4double       truePathLength);
  virtual G4bool IsSingular() const {return true;}
};

#endif
