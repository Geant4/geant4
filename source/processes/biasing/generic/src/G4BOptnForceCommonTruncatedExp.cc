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
#include "G4BOptnForceCommonTruncatedExp.hh"
#include "G4ILawCommonTruncatedExp.hh"
#include "G4ILawForceFreeFlight.hh"
#include "G4TransportationManager.hh"

#include "Randomize.hh"
#include "G4BiasingProcessInterface.hh"

G4BOptnForceCommonTruncatedExp::G4BOptnForceCommonTruncatedExp(G4String name)
  : G4VBiasingOperation(name),
    fNumberOfSharing(0),
    fProcessToApply(nullptr),
    fInteractionOccured(false),
    fMaximumDistance(-1.0)
{
  fCommonTruncatedExpLaw = new G4ILawCommonTruncatedExp("ExpLawForOperation"+name);
  fForceFreeFlightLaw    = new G4ILawForceFreeFlight   ("FFFLawForOperation"+name);
  
  fTotalCrossSection = 0.0;
}

G4BOptnForceCommonTruncatedExp::~G4BOptnForceCommonTruncatedExp()
{
  if ( fCommonTruncatedExpLaw ) delete fCommonTruncatedExpLaw;
  if ( fForceFreeFlightLaw )    delete fForceFreeFlightLaw;
}

const G4VBiasingInteractionLaw* G4BOptnForceCommonTruncatedExp::
ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*       callingProcess, 
				       G4ForceCondition&                proposeForceCondition )
{
  if ( callingProcess->GetWrappedProcess() == fProcessToApply )
    {
      proposeForceCondition =                 Forced;
      return                  fCommonTruncatedExpLaw;
    }
  else
    {
      proposeForceCondition =                 Forced;
      return                     fForceFreeFlightLaw;
    }
}


G4GPILSelection   G4BOptnForceCommonTruncatedExp::ProposeGPILSelection( const G4GPILSelection )
{
  return NotCandidateForSelection;
}


G4VParticleChange* G4BOptnForceCommonTruncatedExp::ApplyFinalStateBiasing( const G4BiasingProcessInterface* callingProcess,
									   const G4Track*                   track,
									   const G4Step*                    step,
									   G4bool&                          forceFinalState )
{
  if ( callingProcess->GetWrappedProcess() != fProcessToApply )
    {
      forceFinalState = true;
      fDummyParticleChange.Initialize( *track );
      return &fDummyParticleChange; 
    }
  if ( fInteractionOccured )
    {
      forceFinalState = true;
      fDummyParticleChange.Initialize( *track );
      return &fDummyParticleChange;
    }
  
  // -- checks if process won the GPIL race:
  G4double processGPIL = callingProcess->GetPostStepGPIL() < callingProcess->GetAlongStepGPIL() ?
    callingProcess->GetPostStepGPIL() : callingProcess->GetAlongStepGPIL() ;
  if ( processGPIL <= step->GetStepLength() )
    {
      // -- if process won, wrapped process produces the final state.
      // -- In this case, the weight for occurence biasing is applied
      // -- by the callingProcess, at exit of present method. This is
      // -- selected by "forceFinalState = false":
      forceFinalState     = false;
      fInteractionOccured = true;
      return callingProcess->GetWrappedProcess()->PostStepDoIt( *track, *step );
    }
  else
    {
      forceFinalState = true;
      fDummyParticleChange.Initialize( *track );
      return &fDummyParticleChange; 
    }
}


void G4BOptnForceCommonTruncatedExp::AddCrossSection( const G4VProcess* process, G4double crossSection )
{
  fTotalCrossSection      += crossSection;
  fCrossSections[process]  = crossSection;
  fNumberOfSharing         = fCrossSections.size();
}


void G4BOptnForceCommonTruncatedExp::Initialize( const G4Track* track )
{
  fCrossSections.clear();
  fTotalCrossSection  = 0.0;
  fNumberOfSharing    = 0;
  fProcessToApply     = 0;
  fInteractionOccured = false;
  fInitialMomentum    = track->GetMomentum();

  G4VSolid* currentSolid = track->GetVolume()->GetLogicalVolume()->GetSolid();
  G4ThreeVector  localPosition = (G4TransportationManager::GetTransportationManager()->
				  GetNavigatorForTracking()->
				  GetGlobalToLocalTransform()).TransformPoint(track->GetPosition());
  G4ThreeVector localDirection = (G4TransportationManager::GetTransportationManager()->
				  GetNavigatorForTracking()->
				  GetGlobalToLocalTransform()).TransformAxis(track->GetMomentumDirection());
  fMaximumDistance = currentSolid->DistanceToOut(localPosition, localDirection);
  if ( fMaximumDistance <= DBL_MIN )  fMaximumDistance = 0.0;
  fCommonTruncatedExpLaw->SetMaximumDistance( fMaximumDistance );
}


void G4BOptnForceCommonTruncatedExp::UpdateForStep( const G4Step* step )
{
  fCrossSections.clear();
  fTotalCrossSection  = 0.0;
  fNumberOfSharing    = 0;
  fProcessToApply     = 0;
  
  fCommonTruncatedExpLaw->UpdateForStep( step->GetStepLength() );
  fMaximumDistance = fCommonTruncatedExpLaw->GetMaximumDistance();
}


void G4BOptnForceCommonTruncatedExp::Sample()
{
  fCommonTruncatedExpLaw->SetForceCrossSection( fTotalCrossSection );
  fCommonTruncatedExpLaw->Sample();
  ChooseProcessToApply();
  fCommonTruncatedExpLaw->SetSelectedProcessXSfraction(fCrossSections[fProcessToApply] / fTotalCrossSection);
}


void G4BOptnForceCommonTruncatedExp::ChooseProcessToApply()
{
  G4double sigmaRand   = G4UniformRand() * fTotalCrossSection;
  G4double sigmaSelect = 0.0;
  for ( std::map< const G4VProcess*, G4double>::const_iterator it = fCrossSections.begin();
	it != fCrossSections.end();
	it++)
    {
      sigmaSelect += (*it).second;
      if ( sigmaRand <= sigmaSelect )
	{
	  fProcessToApply = (*it).first;
	  break;
	}
    }
}
