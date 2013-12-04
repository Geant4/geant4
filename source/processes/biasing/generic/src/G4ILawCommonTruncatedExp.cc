#include "G4ILawCommonTruncatedExp.hh"
#include "G4Track.hh"

#include "G4BOptnForceCommonTruncatedExp.hh"
#include "G4BiasingProcessInterface.hh"

//#include "G4Exception.hh"

G4ILawCommonTruncatedExp::G4ILawCommonTruncatedExp(G4String name)
  : G4VBiasingInteractionLaw(name),
    fExpInteractionLaw("expLawFor"+name)
    //fIsSingular(false)
{}

G4ILawCommonTruncatedExp::~G4ILawCommonTruncatedExp()
{}


G4double G4ILawCommonTruncatedExp::ComputeEffectiveCrossSectionAt(G4double distance) const
{
  return fExpInteractionLaw.ComputeEffectiveCrossSectionAt( distance ) * fOperation->GetTriggeredProcessXSfraction();
}

G4double G4ILawCommonTruncatedExp::ComputeNonInteractionProbabilityAt(G4double distance) const
{
  G4double niProba = fExpInteractionLaw.ComputeNonInteractionProbabilityAt( distance );

  if ( niProba > 0.0 )
    {
      return std::exp( std::log(niProba) / fNumberOfSharing );
    }
  else
    {
      G4ExceptionDescription ed;
      ed << " Negative probability for `" << GetName() << "' p = " << niProba << " distance = " << distance <<  " !!! " << G4endl;
      G4Exception(" G4ILawCommonTruncatedExp::ComputeNonInteractionProbabilityAt(...)",
		  "BIAS.GEN.08",
		  JustWarning,
		  ed);
      return niProba;
    }
}

G4double G4ILawCommonTruncatedExp::SampleInteractionLength()
{
  if ( fFirstSamplingCall ) fInteractionDistance = fExpInteractionLaw.SampleInteractionLength();
  fFirstSamplingCall = false;
  return DBL_MAX;
}

G4double G4ILawCommonTruncatedExp::UpdateInteractionLengthForStep(G4double truePathLength)
{
  if ( fFirstUpdateCall ) 
    fInteractionDistance = fExpInteractionLaw.UpdateInteractionLengthForStep(truePathLength);
  fFirstUpdateCall = false;
  return  DBL_MAX;
}

void G4ILawCommonTruncatedExp::reset()
{
  fFirstSamplingCall       = true;
  fFirstUpdateCall         = true;
}
