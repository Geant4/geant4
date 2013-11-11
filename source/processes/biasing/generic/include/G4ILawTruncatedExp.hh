#ifndef G4ILawTruncatedExp_hh
#define G4ILawTruncatedExp_hh 1

#include "G4VBiasingInteractionLaw.hh"

class G4ILawTruncatedExp : public G4VBiasingInteractionLaw
{
public:
  G4ILawTruncatedExp(G4String name = "expForceInteractionLaw");
  virtual ~G4ILawTruncatedExp();
  
public:
  virtual G4double     ComputeEffectiveCrossSectionAt(G4double               length) const;
  virtual G4double ComputeNonInteractionProbabilityAt(G4double               length) const;
  // -- sample the distribution
  virtual  G4double           SampleInteractionLength();
  // -- move by true path length, this position becomes the new initial point
  virtual G4double     UpdateInteractionLengthForStep(G4double       truePathLength);
  virtual G4bool IsSingular() const {return fIsSingular;}

public:
  void SetForceCrossSection(G4double xs);

public:
  void         SetMaximumDistance(G4double d) { fMaximumDistance = d;}
  G4double     GetMaximumDistance() const     { return fMaximumDistance;}
  G4double GetInteractionDistance() const     { return fInteractionDistance; }

private:
  G4double     fMaximumDistance;
  G4double        fCrossSection;
  G4double fCrossSectionDefined;
  G4bool            fIsSingular;
  G4double fInteractionDistance;

};

#endif
