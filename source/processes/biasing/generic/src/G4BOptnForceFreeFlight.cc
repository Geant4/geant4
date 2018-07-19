//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "G4BOptnForceFreeFlight.hh"
#include "G4ILawForceFreeFlight.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4Step.hh"



G4BOptnForceFreeFlight::G4BOptnForceFreeFlight(G4String name)
  : G4VBiasingOperation    ( name ),
    fCumulatedWeightChange ( -1.0 ),
    fInitialTrackWeight    ( -1.0 ),
    fOperationComplete     ( true )
{
  fForceFreeFlightInteractionLaw = new G4ILawForceFreeFlight("LawForOperation"+name);
}

G4BOptnForceFreeFlight::~G4BOptnForceFreeFlight()
{
  if ( fForceFreeFlightInteractionLaw ) delete fForceFreeFlightInteractionLaw;
}

const G4VBiasingInteractionLaw* G4BOptnForceFreeFlight::ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*, G4ForceCondition& proposeForceCondition )
{
  fOperationComplete = false;
  proposeForceCondition = Forced;
  return fForceFreeFlightInteractionLaw;
}


G4VParticleChange* G4BOptnForceFreeFlight::ApplyFinalStateBiasing( const G4BiasingProcessInterface* callingProcess,
								   const G4Track* track,
								   const G4Step* step,
								   G4bool& forceFinalState)
{
  // -- If the track is reaching the volume boundary, its free flight ends. In this case, its zero
  // -- weight is brought back to non-zero value: its initial weight is restored by the first
  // -- ApplyFinalStateBiasing operation called, and the weight for force free flight is applied
  // -- is applied by each operation.
  // -- If the track is not reaching the volume boundary, it zero weight flight continues.

  fParticleChange.Initialize( *track );
  forceFinalState    = true;
  if ( step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary )
    {
      // -- Sanity checks:
      if ( fInitialTrackWeight <= DBL_MIN )
	{
	  G4ExceptionDescription ed;
	  ed << " Initial track weight is null ! " << G4endl;
	  G4Exception(" G4BOptnForceFreeFlight::ApplyFinalStateBiasing(...)",
		      "BIAS.GEN.05",
		      JustWarning,
		      ed);
	}
      if ( fCumulatedWeightChange <= DBL_MIN )
	{
	  G4ExceptionDescription ed;
	  ed << " Cumulated weight is null ! " << G4endl;
	  G4Exception(" G4BOptnForceFreeFlight::ApplyFinalStateBiasing(...)",
		      "BIAS.GEN.06",
		      JustWarning,
		      ed);
	}

      G4double proposedWeight = track->GetWeight();
      if ( callingProcess->GetIsFirstPostStepDoItInterface() ) proposedWeight  = fCumulatedWeightChange * fInitialTrackWeight;
      else                                                     proposedWeight *= fCumulatedWeightChange;
      fParticleChange.ProposeWeight(proposedWeight);
      fOperationComplete = true;
    }
  
  return &fParticleChange;
}


void G4BOptnForceFreeFlight::AlongMoveBy( const G4BiasingProcessInterface*, const G4Step*, G4double weightChange )
{
  fCumulatedWeightChange *= weightChange;
}
