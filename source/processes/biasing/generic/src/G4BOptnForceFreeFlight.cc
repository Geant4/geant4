#include "G4BOptnForceFreeFlight.hh"
#include "G4ILawForceFreeFlight.hh"
#include "G4Step.hh"



G4BOptnForceFreeFlight::G4BOptnForceFreeFlight(G4String name)
  : G4VBiasingOperation(name)
{
  fForceFreeFlightInteractionLaw = new G4ILawForceFreeFlight("LawForOperation"+name);
}

G4BOptnForceFreeFlight::~G4BOptnForceFreeFlight()
{}

const G4VBiasingInteractionLaw* G4BOptnForceFreeFlight::ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface* )
{
  return fForceFreeFlightInteractionLaw;
}

G4bool G4BOptnForceFreeFlight::DenyProcessPostStepDoIt( const G4BiasingProcessInterface*, const G4Track*, const G4Step* step, G4double& proposedWeight )
{
  // -- force free flight always deny process to apply its doit.
  // -- if reaching boundary, track is restored with non-zero weight
  if ( fInitialTrackWeight <= DBL_MIN )
    {
      G4ExceptionDescription ed;
      ed << " Initial track weight is null ! " << G4endl;
      G4Exception(" G4BOptnForceFreeFlight::DenyProcessPostStepDoIt(...)",
		  "BIAS.GEN.05",
		  JustWarning,
		  ed);
    }
  if ( fCumulatedWeightChange <= DBL_MIN )
    {
      G4ExceptionDescription ed;
      ed << " Cumulated weight is null ! " << G4endl;
      G4Exception(" G4BOptnForceFreeFlight::DenyProcessPostStepDoIt(...)",
		  "BIAS.GEN.06",
		  JustWarning,
		  ed);
    }
  if ( step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary )
    {
      if ( proposedWeight <= DBL_MIN ) proposedWeight  = fCumulatedWeightChange * fInitialTrackWeight;
      else                             proposedWeight *= fCumulatedWeightChange;
    }
  
  return true;
}

void G4BOptnForceFreeFlight::AlongMoveBy( const G4BiasingProcessInterface*, const G4Step*, G4double weightChange )
{
  fCumulatedWeightChange *= weightChange;
}
