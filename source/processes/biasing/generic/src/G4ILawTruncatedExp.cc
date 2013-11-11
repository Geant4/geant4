#include "G4ILawTruncatedExp.hh"
#include "Randomize.hh"
#include "G4Track.hh"

//#include "G4Exception.hh"

G4ILawTruncatedExp::G4ILawTruncatedExp(G4String name)
  : G4VBiasingInteractionLaw(name),
    fMaximumDistance(0.0),
    fCrossSection(0.0),
    fCrossSectionDefined(false),
    fIsSingular(false)
{}

G4ILawTruncatedExp::~G4ILawTruncatedExp()
{}

void G4ILawTruncatedExp::SetForceCrossSection(G4double crossSection)
{
  if (crossSection < 0.0)
    {
      G4Exception("G4ILawTruncatedExp::SetForceCrossSection(..)",
		  "NEGATIVE_XS",
		  JustWarning,
		  "Cross-section value passed is negative. It is set to zero.");
      fIsSingular  = true;
      crossSection = 0.0;
    }
  fIsSingular          = false;
  fCrossSectionDefined = true;
  fCrossSection        = crossSection;
}

G4double G4ILawTruncatedExp::ComputeEffectiveCrossSectionAt(G4double distance) const
{
  if ( !fCrossSectionDefined )
    {
      G4Exception("G4ILawTruncatedExp::ComputeEffectiveCrossSection(..)",
		  "UNDEFINED_XS",
		  JustWarning,
		  "Cross-section value requested, but has not been defined yet. Assumes 0.");
      // -- zero cross-section, returns the limit form of the effective cross-section:
      return 1.0 / ( fMaximumDistance - distance );
    }
  G4double denum = 1.0 - exp( -fCrossSection * ( fMaximumDistance - distance) );
  return fCrossSection / denum;
}

G4double G4ILawTruncatedExp::ComputeNonInteractionProbabilityAt(G4double distance) const
{
  if (!fCrossSectionDefined)
    {
      G4Exception("G4ILawTruncatedExp::ComputeNonInteractionProbability(..)",
		  "UNDEFINED_XS",
		  JustWarning,
		  "Non interaction probabitlity value requested, but cross section has not been defined yet. Assumes it to be 0.");
      // -- return limit case of null cross-section:
      return 1.0 - distance / fMaximumDistance;
    }
  G4double   num = 1.0 - exp( -fCrossSection*distance);
  G4double denum = 1.0 - exp( -fCrossSection*fMaximumDistance);
  return 1.0 - num/denum;
}

G4double G4ILawTruncatedExp::SampleInteractionLength()
{
  if ( !fCrossSectionDefined )
    {
      G4Exception("G4ILawTruncatedExp::Sample(..)",
		  "UNDEFINED_XS",
		  JustWarning,
		  "Trying to sample while cross-section is not defined, assuming 0.");
      fInteractionDistance = G4UniformRand() * fMaximumDistance;
      return fInteractionDistance;
    }
  fInteractionDistance = -log(1.0 - G4UniformRand()* (1.0 - exp(-fCrossSection*fMaximumDistance)))/fCrossSection;
  return fInteractionDistance;
}


G4double G4ILawTruncatedExp::UpdateInteractionLengthForStep(G4double       truePathLength)
{
  fInteractionDistance -= truePathLength;
  fMaximumDistance     -= truePathLength;

  if ( fInteractionDistance < 0 )  G4cout << GetName() << " 1 " <<
				     " NEGATIVE # of int length !!!!!!!!!!!!!!! MESSAGE TO BE CHANGED !!!!!!!!!!!!! " 
					  << fInteractionDistance << G4endl;
  return  fInteractionDistance;
}
