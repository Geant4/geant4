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
  if ( fInitialTrackWeight <= DBL_MIN ) G4cout << " 1. Houston, we get a problem : initial weight ZERO " << G4endl;
  if ( fCumulatedWeightChange <= DBL_MIN ) G4cout << " 2. Houston, we get a problem : cumulated weight ZERO " << G4endl;
  if ( step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary )
    {
      if ( proposedWeight <= DBL_MIN ) proposedWeight  = fCumulatedWeightChange * fInitialTrackWeight;
      else                             proposedWeight *= fCumulatedWeightChange;
    }
  //  G4cout << " ----- DenyProcessPostStepDoIt " << fCumulatedWeightChange << " "  << proposedWeight << " " <<  1. - proposedWeight<< G4endl; // !
  return true;
}

void G4BOptnForceFreeFlight::AlongMoveBy( const G4BiasingProcessInterface*, const G4Step*, G4double weightChange )
{
  //  G4cout << " ----- fCumulatedWeightChange , weightChange " <<  fCumulatedWeightChange << " " << weightChange << " " << this << G4endl; // !
  fCumulatedWeightChange *= weightChange;
}
