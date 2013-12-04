#include "G4InteractionLawPhysical.hh"
#include "Randomize.hh"

G4InteractionLawPhysical::G4InteractionLawPhysical(G4String name)
  : G4VBiasingInteractionLaw(name),
    fCrossSection(0.0),
    fCrossSectionDefined(false),
    fNumberOfInteractionLength(-1.0)
{}

G4InteractionLawPhysical::~G4InteractionLawPhysical()
{}

void G4InteractionLawPhysical::SetPhysicalCrossSection(G4double crossSection)
{
  if (crossSection < 0.0)
    {
      G4Exception("G4InteractionLawPhysical::SetPhysicalCrossSection(..)",
		  "BIAS.GEN.14",
		  JustWarning,
		  "Cross-section value passed is negative. It is set to zero !");
      crossSection = 0.0;
    }
  fCrossSectionDefined = true;
  fCrossSection        = crossSection;
}

G4double G4InteractionLawPhysical::ComputeEffectiveCrossSectionAt(G4double) const
{
  if (!fCrossSectionDefined) G4Exception("G4InteractionLawPhysical::ComputeEffectiveCrossSection(..)",
					 "BIAS.GEN.15",
					 JustWarning,
					 "Cross-section value requested, but has not been defined yet. Assumes 0 !");
  return fCrossSection;
}

G4double G4InteractionLawPhysical::ComputeNonInteractionProbabilityAt(G4double stepLength) const
{
  if (!fCrossSectionDefined) G4Exception("G4InteractionLawPhysical::ComputeNonInteractionProbability(..)",
					 "BIAS.GEN.16",
					 JustWarning,
					 "Non interaction probabitlity value requested, but cross section has not been defined yet. Assumes it to be 0 !");
  // -- allows zero cross-section case, by convention:
  if ( fCrossSection == 0.0 ) return 1.0;
  else return std::exp(-fCrossSection*stepLength);
}

G4double G4InteractionLawPhysical::SampleInteractionLength()
{
  if ( !fCrossSectionDefined || fCrossSection < 0.0 )  G4Exception("G4InteractionLawPhysical::Sample(..)",
								   "BIAS.GEN.17",
								   FatalException,
								   "Trying to sample while cross-section is not defined or < 0 !");
  if ( fCrossSection == 0.0 ) return DBL_MAX;

  fNumberOfInteractionLength =  -std::log( G4UniformRand() );
  return fNumberOfInteractionLength/fCrossSection;
}


G4double G4InteractionLawPhysical::UpdateInteractionLengthForStep(G4double       truePathLength)
{
  fNumberOfInteractionLength -= truePathLength*fCrossSection;
  if ( fNumberOfInteractionLength < 0 ) 
    {
      G4ExceptionDescription ed;
      ed << " Negative number of interaction length for `" << GetName() << "' " << fNumberOfInteractionLength << ", set it to zero !" << G4endl; 
      G4Exception("G4InteractionLawPhysical::UpdateInteractionLengthForStep(...)",
		  "BIAS.GEN.13",
		  JustWarning,
		  ed);
      fNumberOfInteractionLength = 0.;
    }
  return  fNumberOfInteractionLength;
}
