#include "G4ILawForceFreeFlight.hh"

G4ILawForceFreeFlight::G4ILawForceFreeFlight(G4String name)
  : G4VBiasingInteractionLaw(name)
{}

G4ILawForceFreeFlight::~G4ILawForceFreeFlight()
{}

G4double G4ILawForceFreeFlight::ComputeEffectiveCrossSectionAt(G4double) const
{
  return 0.0;
}

G4double G4ILawForceFreeFlight::ComputeNonInteractionProbabilityAt(G4double) const
{
  return 1.0;
}

G4double G4ILawForceFreeFlight::SampleInteractionLength()
{
  return DBL_MAX;
}

G4double G4ILawForceFreeFlight::UpdateInteractionLengthForStep(G4double)
{
  return DBL_MAX;
}
