#ifndef G4InteractionLawPhysical_hh
#define G4InteractionLawPhysical_hh 1

#include "G4VBiasingInteractionLaw.hh"

class G4InteractionLawPhysical : public G4VBiasingInteractionLaw
{
public:
  G4InteractionLawPhysical(G4String name = "exponentialLaw");
  virtual ~G4InteractionLawPhysical();

public:
  void     SetPhysicalCrossSection(G4double crossSection);
  G4double GetPhysicalCrossSection() const {return fCrossSection;}
  
public:
  virtual void                     DefineInitialPoint(const G4Track*         track);
  virtual G4double     ComputeEffectiveCrossSectionAt(G4double               length) const;
  virtual G4double ComputeNonInteractionProbabilityAt(G4double               length) const;
  // -- sample the distribution
  virtual  G4double           SampleInteractionLength();
  // -- move by true path length, this position becomes the new initial point
  virtual G4double     UpdateInteractionLengthForStep(G4double       truePathLength);


private:
  G4double fCrossSection;
  G4bool   fCrossSectionDefined;
  G4double fNumberOfInteractionLength;

};

#endif
