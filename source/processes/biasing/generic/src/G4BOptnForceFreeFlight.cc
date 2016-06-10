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
