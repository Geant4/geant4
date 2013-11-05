#include "G4InteractionLawPhysical.hh"
#include "Randomize.hh"

//#include "G4Exception.hh"

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
		  "NEGATIVE_XS",
		  JustWarning,
		  "Cross-section value passed is negative. It is set to zero.");
      crossSection = 0.0;
    }
  fCrossSectionDefined = true;
  fCrossSection        = crossSection;
}

G4double G4InteractionLawPhysical::ComputeEffectiveCrossSectionAt(G4double) const
{
  if (!fCrossSectionDefined) G4Exception("G4InteractionLawPhysical::ComputeEffectiveCrossSection(..)",
					 "UNDEFINED_XS",
					 JustWarning,
					 "Cross-section value requested, but has not been defined yet. Assumes 0.");
  return fCrossSection;
}

G4double G4InteractionLawPhysical::ComputeNonInteractionProbabilityAt(G4double stepLength) const
{
  if (!fCrossSectionDefined) G4Exception("G4InteractionLawPhysical::ComputeNonInteractionProbability(..)",
					 "UNDEFINED_XS",
					 JustWarning,
					 "Non interaction probabitlity value requested, but cross section has not been defined yet. Assumes it to be 0.");
  // -- allows zero cross-section case, by convention:
  if ( fCrossSection == 0.0 ) return 1.0;
  else return exp(-fCrossSection*stepLength);
}

G4double G4InteractionLawPhysical::SampleInteractionLength()
{
  if ( !fCrossSectionDefined || fCrossSection < 0.0 )  G4Exception("G4InteractionLawPhysical::Sample(..)",
								   "UNDEFINED_XS",
								   FatalException,
								   "Trying to sample while cross-section is not defined or < 0.");
  if ( fCrossSection == 0.0 ) return DBL_MAX;

  fNumberOfInteractionLength =  -std::log( G4UniformRand() );
  return fNumberOfInteractionLength/fCrossSection;
}


G4double G4InteractionLawPhysical::UpdateInteractionLengthForStep(G4double       truePathLength)
{
  fNumberOfInteractionLength -= truePathLength*fCrossSection;
  if ( fNumberOfInteractionLength < 0 )  G4cout << GetName() <<
					   " NEGATIVE # of int length !!!!!!!!!!!!!!! MESSAGE TO BE CHANGED !!!!!!!!!!!!! " 
						<< fNumberOfInteractionLength << G4endl;
  return  fNumberOfInteractionLength;
}


void G4InteractionLawPhysical::DefineInitialPoint(const G4Track*)
{}
