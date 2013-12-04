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
#include "G4BiasingProcessInterface.hh"
#include "G4VBiasingOperator.hh"
#include "G4VBiasingOperation.hh"
#include "G4ParticleChangeForOccurenceBiasing.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleChangeForNothing.hh"
#include "G4VBiasingInteractionLaw.hh"
#include "G4InteractionLawPhysical.hh"
#include "G4ProcessManager.hh"
#include "G4BiasingTrackData.hh"
#include "G4BiasingTrackDataStore.hh"
#include "G4BiasingAppliedCase.hh"

G4Cache<G4bool>                     G4BiasingProcessInterface::fResetInteractionLaws;// = true;
G4Cache<G4bool>                     G4BiasingProcessInterface::fCommonStart;//          = true;
G4Cache<G4bool>                     G4BiasingProcessInterface::fCommonEnd;//            = true;
G4MapCache< const G4ProcessManager*, std::vector< G4BiasingProcessInterface* > > G4BiasingProcessInterface::fManagerInterfaceMap;


G4BiasingProcessInterface::G4BiasingProcessInterface(G4String name)
  :  G4VProcess( name ),
     fCurrentBiasingOperator ( 0 ),
     fPreviousBiasingOperator( 0 ),
     fWrappedProcess        ( 0     ),
     fIsPhysicsBasedBiasing ( false ),
     fWrappedProcessIsAtRest( false ),
     fWrappedProcessIsAlong ( false ),
     fWrappedProcessIsPost  ( false ),
     fWrappedProcessInteractionLength( -1.0 ),
     fBiasingInteractionLaw ( 0 ),
     fPhysicalInteractionLaw( 0 ),
     fOccurenceBiasingParticleChange( 0 ),
     fIamFirstGPIL          ( false )
{
  for (G4int i = 0 ; i < 8 ; i++)  fFirstLastFlags[i] = false;
    fResetInteractionLaws.Put( true );
    fCommonStart.Put(true);
    fCommonEnd.Put(true);
}


G4BiasingProcessInterface::G4BiasingProcessInterface(G4VProcess* wrappedProcess,
						     G4bool wrappedIsAtRest, G4bool wrappedIsAlongStep, G4bool wrappedIsPostStep,
						     G4String useThisName)
  : G4VProcess( useThisName != "" ? useThisName : "biasWrapper("+wrappedProcess->GetProcessName()+")",
	       wrappedProcess->GetProcessType()),
    fCurrentBiasingOperator ( 0 ),
    fPreviousBiasingOperator( 0 ),
    fWrappedProcess        ( wrappedProcess     ),
    fIsPhysicsBasedBiasing ( true               ),
    fWrappedProcessIsAtRest( wrappedIsAtRest    ),
    fWrappedProcessIsAlong ( wrappedIsAlongStep ),
    fWrappedProcessIsPost  ( wrappedIsPostStep  ),
    fWrappedProcessInteractionLength( -1.0 ),
    fBiasingInteractionLaw ( 0 ),
    fPhysicalInteractionLaw( 0 ),
    fOccurenceBiasingParticleChange( 0 ),
    fIamFirstGPIL          ( false )
{
  for (G4int i = 0 ; i < 8 ; i++)  fFirstLastFlags[i] = false;
  
  SetProcessSubType(fWrappedProcess->GetProcessSubType());

  // -- create physical interaction law:
  fPhysicalInteractionLaw         = new            G4InteractionLawPhysical("PhysicalInteractionLawFor("+GetProcessName()+")");
  // -- instantiate particle change wrapper for occurence biaising:
  fOccurenceBiasingParticleChange = new G4ParticleChangeForOccurenceBiasing("biasingPCfor"+GetProcessName());
  fParticleChange                 = new                    G4ParticleChange();
  // -- instantiate a "do nothing" particle change:
  fDummyParticleChange            = new          G4ParticleChangeForNothing();
}



G4BiasingProcessInterface::~G4BiasingProcessInterface()
{
  if ( fPhysicalInteractionLaw  != 0   ) delete fPhysicalInteractionLaw;
  if ( fOccurenceBiasingParticleChange ) delete fOccurenceBiasingParticleChange;
  if ( fDummyParticleChange            ) delete fDummyParticleChange;
}


void G4BiasingProcessInterface::StartTracking(G4Track* track)
{
  fCurrentTrack = track;
  if ( fIsPhysicsBasedBiasing ) fWrappedProcess->StartTracking(fCurrentTrack);
  fCurrentBiasingOperator             = 0;
  fPreviousBiasingOperator            = 0;
  fOccurenceBiasingOperation          = 0;
  fPreviousOccurenceBiasingOperation  = 0;
  fFinalStateBiasingOperation         = 0;
  fPreviousFinalStateBiasingOperation = 0;
  fNonPhysicsBiasingOperation         = 0;
  fPreviousNonPhysicsBiasingOperation = 0;
  fBiasingInteractionLaw              = 0;
  fPreviousBiasingInteractionLaw      = 0;
  
  fPreviousStepSize           = -1.0;
  
  
  fResetWrappedProcessInteractionLength = false;
  
  if ( fCommonStart.Get() )
    {
        fCommonStart.Put( false );// = false;
        fCommonEnd.Put(true);//   = true;
      
      for ( size_t optr = 0 ; optr < ( G4VBiasingOperator::GetBiasingOperators() ).size() ; optr ++)
	( G4VBiasingOperator::GetBiasingOperators() )[optr]->StartTracking( fCurrentTrack );
    }
}


void G4BiasingProcessInterface::EndTracking()
{
  if ( fIsPhysicsBasedBiasing ) fWrappedProcess->EndTracking();
  if ( fCurrentBiasingOperator) fCurrentBiasingOperator->ExitingBiasing( fCurrentTrack, this ); 
  fCurrentBiasingOperator   = 0;
  fPreviousBiasingOperator  = 0;
  fBiasingInteractionLaw = 0;
  

  // -- !! this part might have to be improved : could be time consuming
  // -- !! and assumes all tracks are killed during tracking, which is
  // -- !! not true : stacking operations may kill tracks
  if ( fCommonEnd.Get() )
    {
        fCommonEnd.Put( false );//   = false;
        fCommonStart.Put( true );//  = true;
      
      for ( size_t optr = 0 ; optr < ( G4VBiasingOperator::GetBiasingOperators() ).size() ; optr ++)
	( G4VBiasingOperator::GetBiasingOperators() )[optr]->EndTracking( );
      
      if ( ( fCurrentTrack->GetTrackStatus() == fStopAndKill ) || ( fCurrentTrack->GetTrackStatus() == fKillTrackAndSecondaries ) )
	{
	  G4BiasingTrackData* biasingData = G4BiasingTrackDataStore::GetInstance()->GetBiasingTrackData( fCurrentTrack );
	  if ( biasingData ) delete biasingData; // -- this also deregisters the biasing data from the track data store
	  if ( fCurrentTrack->GetTrackStatus() == fKillTrackAndSecondaries )
	    {
	      const G4TrackVector* secondaries = fCurrentTrack->GetStep()->GetSecondary();
	      for ( size_t i2nd = 0 ; i2nd < secondaries->size() ; i2nd++ )
		{
		  biasingData = G4BiasingTrackDataStore::GetInstance()->GetBiasingTrackData( (*secondaries)[i2nd] );
		  if ( biasingData ) delete biasingData;
		}
	    }
	}
    }
}



G4double G4BiasingProcessInterface::PostStepGetPhysicalInteractionLength(const G4Track&                track,
									 G4double           previousStepSize,
									 G4ForceCondition*         condition)
{
  // -- Remember previous operator and proposed operations, if any, and reset:
  // -------------------------------------------------------------------------
  // -- remember:
  fPreviousBiasingOperator            =     fCurrentBiasingOperator;
  fPreviousOccurenceBiasingOperation  =  fOccurenceBiasingOperation;
  fPreviousFinalStateBiasingOperation = fFinalStateBiasingOperation;
  fPreviousNonPhysicsBiasingOperation = fNonPhysicsBiasingOperation;
  fPreviousBiasingInteractionLaw      =      fBiasingInteractionLaw;
  // -- reset:
  fOccurenceBiasingOperation          =                           0;
  fFinalStateBiasingOperation         =                           0;
  fNonPhysicsBiasingOperation         =                           0;
  fBiasingInteractionLaw              =                           0;
  // -- Physics PostStep and AlongStep GPIL
  fWrappedProcessPostStepGPIL         =                     DBL_MAX;
  fBiasingPostStepGPIL                =                     DBL_MAX;
  fWrappedProcessInteractionLength    =                     DBL_MAX; // -- inverse of analog cross-section, no biasing counter-part in general
  fWrappedProcessForceCondition       =                   NotForced;
  fBiasingForceCondition              =                   NotForced;
  fWrappedProcessAlongStepGPIL        =                     DBL_MAX;
  fBiasingAlongStepGPIL               =                     DBL_MAX;
  fWrappedProcessGPILSelection        =    NotCandidateForSelection;
  fBiasingGPILSelection               =    NotCandidateForSelection;
  // -- for helper:
  fPreviousStepSize                   =            previousStepSize;


  // -- If new volume, get possible new biasing operator:
  // ----------------------------------------------------
  // -- [Note : bug with this first step ! Does not work if previous step was concurrently limited with geometry]
  G4bool  firstStepInVolume = ( (track.GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) || (track.GetCurrentStepNumber() == 1) );
  if ( firstStepInVolume ) fCurrentBiasingOperator = G4VBiasingOperator::GetBiasingOperator( track.GetVolume()->GetLogicalVolume() );

  // -----------------------------------
  // -- Previous step was under biasing:
  // -----------------------------------
  if ( fPreviousBiasingOperator != 0 )
    {
      // -- if current step not under same operator, let this operator knows:
      if ( fPreviousBiasingOperator != fCurrentBiasingOperator ) fPreviousBiasingOperator->ExitingBiasing( &track, this ); 
      // -- if biasing does not continue, set process behavior to standard tracking:
      // -- ... and see [*] below...
      if ( fCurrentBiasingOperator == 0 )
	{
	  fResetWrappedProcessInteractionLength = true;
	  ResetForUnbiasedTracking();
	}
    }
  
  // --------------------------------------------------------------
  // -- no operator : analog tracking if physics-based, or nothing:
  // --------------------------------------------------------------
  if ( fCurrentBiasingOperator == 0 )
    {
      if ( fIsPhysicsBasedBiasing ) 
	{
	  // -- [*] this wrapped process has just been ResetForUnbiasedTracking(), it has a fresh #int length, and is
	  // -- in a state where it believes it is starting tracking : let it believe hence the previous step
	  // -- length was 0.0, as for a normal first step:
	  if ( fResetWrappedProcessInteractionLength )
	    {
	      fResetWrappedProcessInteractionLength = false;
	      fWrappedProcess->ResetNumberOfInteractionLengthLeft();
	      return fWrappedProcess->PostStepGetPhysicalInteractionLength(track, 0.0, condition);
	    }
	  return fWrappedProcess->PostStepGetPhysicalInteractionLength(track, previousStepSize, condition);
	}
      else
  	{
  	  *condition = NotForced;
  	  return DBL_MAX;
  	}
    }

  

  // --------------------------------------------------
  // -- An biasing operator exists. Proceed with
  // -- treating non-physics and physics biasing cases:
  //---------------------------------------------------
  
  // -- non-physics-based biasing case:
  // ----------------------------------
  if ( !fIsPhysicsBasedBiasing )
    {  
      fNonPhysicsBiasingOperation = fCurrentBiasingOperator->GetProposedNonPhysicsBiasingOperation( &track, this );
      if ( fNonPhysicsBiasingOperation == 0 )
	{
	  *condition = NotForced;
	  return DBL_MAX;
	}
      return fNonPhysicsBiasingOperation->DistanceToApplyOperation(&track, previousStepSize, condition);
    }
  

  // -- Physics based biasing case:
  // ------------------------------
  // -- call to underneath physics process PostStepGPIL to update it with current point:
  fWrappedProcessPostStepGPIL   = fWrappedProcess->PostStepGetPhysicalInteractionLength(track, previousStepSize, condition);
  fWrappedProcessForceCondition = *condition;
  // -- **! At this point, might have to be careful of previousStepSize being larger than the previous process
  // -- **! PostStepGPIL proposed: as we disregard this PostStepGPIL value, such situation does happen. Does not
  // -- **! look to have generated problems for now, but might be fragile.

  // -- Ask for possible GPIL biasing operation:
  fOccurenceBiasingOperation = fCurrentBiasingOperator->GetProposedOccurenceBiasingOperation( &track, this );

  // -- no operation for occurence biasing, analog GPIL returns the wrapped process GPIL and condition values
  // -- (note that condition was set above):
  if ( fOccurenceBiasingOperation == 0 ) return fWrappedProcessPostStepGPIL;

  // -- A valid GPIL biasing operation has been proposed:
  // -- 0) remember wrapped process will need to be reset on biasing exit, if particle survives:
  fResetWrappedProcessInteractionLength = true;
  // -- 1) collect/update process interaction length for reference analog interaction law:
  fWrappedProcessInteractionLength = fWrappedProcess->GetCurrentInteractionLength();
  fPhysicalInteractionLaw->SetPhysicalCrossSection( 1.0 / fWrappedProcessInteractionLength );
  // -- 2) Collect biasing interaction law:
  // --    The interaction law pointer is collected as a const pointer to the interaction law object.
  // --    This interaction law will be kept under control of the biasing operation, which is the only
  // --    entity that will change the state of the biasing interaction law.
  fBiasingInteractionLaw            = fOccurenceBiasingOperation->ProvideOccurenceBiasingInteractionLaw( this );
  // -- 3) Ask operation to sample the biasing interaction law:
  fBiasingPostStepGPIL              = fBiasingInteractionLaw->GetSampledInteractionLength();
  fBiasingForceCondition            = fOccurenceBiasingOperation->ProposeForceCondition( fWrappedProcessForceCondition );

  // -- finish
  *condition = fBiasingForceCondition;
  return       fBiasingPostStepGPIL;

}



G4VParticleChange* G4BiasingProcessInterface::PostStepDoIt(const G4Track& track,
							   const G4Step&   step)
{
  // ---------------------------------------
  // -- case outside of volume with biasing:
  // ---------------------------------------
  if ( fCurrentBiasingOperator == 0 ) return fWrappedProcess->PostStepDoIt(track, step);
  
  // ----------------------------
  // -- non-physics biasing case:
  // ----------------------------
  if ( !fIsPhysicsBasedBiasing )
    {
      G4VParticleChange* particleChange = fNonPhysicsBiasingOperation->GenerateBiasingFinalState( &track, &step );
      fCurrentBiasingOperator->ReportOperationApplied( this, BAC_NonPhysics, fNonPhysicsBiasingOperation, particleChange );
      return particleChange;
    }

  // -- physics biasing case:
  // ------------------------
  // -- It proceeds with the following logic:
  // -- 1) If an occurence biasing operation exists, it makes the
  // --    decision about the interaction to happen or not.
  // --    If the occurence operation refuses the interaction, it
  // --    can only propose a new track weight.
  // -- 2) The interaction is produced by the final state biasing
  // --    operation, if it exists, or by the wrapped process in
  // --    the other case.
  // -- Hence 2) happens if 1) decides so
  //          2) happens if there is no occurence biasing operation
  if ( fOccurenceBiasingOperation != 0 )
    {
      G4double proposedTrackWeight = track.GetWeight();
      if ( fOccurenceBiasingOperation->DenyProcessPostStepDoIt( this, &track, &step, proposedTrackWeight ) )
	{
	  fParticleChange->Initialize( track ); // **??** <= might  use a light version for particle change here
	  fParticleChange->ProposeParentWeight( proposedTrackWeight );
	  fCurrentBiasingOperator->ReportOperationApplied( this, BAC_DenyInteraction, fOccurenceBiasingOperation, fParticleChange );
	  return fParticleChange;
	}
    }

  
  // -- usual case with generated final state:
  G4VParticleChange*   finalStateParticleChange;
  G4BiasingAppliedCase BAC;
  fFinalStateBiasingOperation = fCurrentBiasingOperator->GetProposedFinalStateBiasingOperation( &track, this );
  if ( fFinalStateBiasingOperation != 0 )
    {
      finalStateParticleChange = fFinalStateBiasingOperation->ApplyFinalStateBiasing( this, &track, &step );
      BAC = BAC_FinalState;
    }
  else
    {
      finalStateParticleChange = fWrappedProcess->PostStepDoIt(track, step);
      BAC =  BAC_None ;
    }
  // -- if no occurence biasing operation, we're done:
  if ( fOccurenceBiasingOperation == 0 )
    {
      fCurrentBiasingOperator->ReportOperationApplied( this, BAC, fFinalStateBiasingOperation, finalStateParticleChange );
      return finalStateParticleChange;
    }
  
  // -- If occurence biasing, applies on top of final state occurence biasing weight correction:
  G4double weightForInteraction = 1.0;
  if ( !fBiasingInteractionLaw->IsSingular() ) weightForInteraction =
						 fPhysicalInteractionLaw->ComputeEffectiveCrossSectionAt(step.GetStepLength()) /
						 fBiasingInteractionLaw ->ComputeEffectiveCrossSectionAt(step.GetStepLength());
  else
    {
      // -- at this point effective XS can only be infinite, if not, there is a logic problem
      if ( !fBiasingInteractionLaw->IsEffectiveCrossSectionInfinite() )
	{
	  G4ExceptionDescription ed;
	  ed << "Internal inconsistency in cross-section handling. Please report !" << G4endl;
	  G4Exception(" G4BiasingProcessInterface::PostStepDoIt(...)",
		      "BIAS.GEN.02",
		      JustWarning,
		      ed);
	  // -- if XS is infinite, weight is zero (and will stay zero), but we'll do differently.
	  // -- Should foresee in addition something to remember that in case of singular
	  // -- distribution, weight can only be partly calculated
	}
    }
  
  if ( weightForInteraction <= 0. )
    {
      G4ExceptionDescription ed;
      ed << " Negative interaction weight : w_I = "
	 <<  weightForInteraction <<
	" XS_I(phys) = " << fBiasingInteractionLaw ->ComputeEffectiveCrossSectionAt(step.GetStepLength()) <<
	" XS_I(bias) = " << fPhysicalInteractionLaw->ComputeEffectiveCrossSectionAt(step.GetStepLength()) <<
	" step length = " << step.GetStepLength() <<
	" Interaction law = `" << fBiasingInteractionLaw << "'" <<
	G4endl;
      G4Exception(" G4BiasingProcessInterface::PostStepDoIt(...)",
		  "BIAS.GEN.03",
		  JustWarning,
		  ed);
      
    }
  
  fCurrentBiasingOperator->ReportOperationApplied( this,                        BAC,
						   fOccurenceBiasingOperation,  weightForInteraction,
						   fFinalStateBiasingOperation, finalStateParticleChange );  

  fOccurenceBiasingParticleChange->SetOccurenceWeightForInteraction( weightForInteraction );
  fOccurenceBiasingParticleChange->SetSecondaryWeightByProcess( true );
  fOccurenceBiasingParticleChange->SetWrappedParticleChange( finalStateParticleChange );
  fOccurenceBiasingParticleChange->ProposeTrackStatus( finalStateParticleChange->GetTrackStatus() );
  fOccurenceBiasingParticleChange->StealSecondaries(); // -- this also makes weightForInteraction applied to secondaries stolen

  // -- finish:
  return fOccurenceBiasingParticleChange;
  
}


// -- AlongStep methods:
G4double           G4BiasingProcessInterface::AlongStepGetPhysicalInteractionLength(const G4Track&                track,
										    G4double           previousStepSize,
										    G4double         currentMinimumStep, 
										    G4double&            proposedSafety, 
										    G4GPILSelection*          selection)
{
  // -- for helper methods:
  fCurrentMinimumStep = currentMinimumStep;
  fProposedSafety     = proposedSafety;
  

  // -- initialization default case:
  fWrappedProcessAlongStepGPIL =                  DBL_MAX;
  *selection                   = NotCandidateForSelection;
  // ---------------------------------------
  // -- case outside of volume with biasing:
  // ---------------------------------------
  if ( fCurrentBiasingOperator == 0 )
    {
      if ( fWrappedProcessIsAlong ) fWrappedProcessAlongStepGPIL = 
				      fWrappedProcess->AlongStepGetPhysicalInteractionLength(track,
											     previousStepSize,
											     currentMinimumStep, 
											     proposedSafety, 
											     selection);
      return fWrappedProcessAlongStepGPIL;
    }
  
  // --------------------------------------------------------------------
  // -- non-physics based biasing: no along operation expected (for now):
  // --------------------------------------------------------------------
  if ( !fIsPhysicsBasedBiasing ) return fWrappedProcessAlongStepGPIL;
  
  // ----------------------
  // -- physics-based case:
  // ----------------------
  if ( fOccurenceBiasingOperation == 0 )
    {
      if ( fWrappedProcessIsAlong ) fWrappedProcessAlongStepGPIL =
				      fWrappedProcess->AlongStepGetPhysicalInteractionLength(track,
											     previousStepSize,
											     currentMinimumStep, 
											     proposedSafety, 
											     selection);
      return fWrappedProcessAlongStepGPIL;
    }

  
  // ----------------------------------------------------------
  // -- From here we have an valid occurence biasing operation:
  // ----------------------------------------------------------
  // -- Give operation opportunity to shorten step proposed by physics process:
  fBiasingAlongStepGPIL =  fOccurenceBiasingOperation->ProposeAlongStepLimit( this );
  G4double minimumStep  = fBiasingAlongStepGPIL < currentMinimumStep ? fBiasingAlongStepGPIL : currentMinimumStep ;
  // -- wrapped process is called with minimum step ( <= currentMinimumStep passed ) : an along process can not
  // -- have its operation stretched over what it expects:
  if ( fWrappedProcessIsAlong )
    {
      fWrappedProcessAlongStepGPIL = fWrappedProcess->AlongStepGetPhysicalInteractionLength(track,
											    previousStepSize,
											    minimumStep, 
											    proposedSafety, 
											    selection);
      fWrappedProcessGPILSelection = *selection;
      fBiasingGPILSelection        = fOccurenceBiasingOperation->ProposeGPILSelection( fWrappedProcessGPILSelection );
    }
  else
    {
      fBiasingGPILSelection        = fOccurenceBiasingOperation->ProposeGPILSelection( NotCandidateForSelection );
      fWrappedProcessAlongStepGPIL = fBiasingAlongStepGPIL;
    }
  
  *selection = fBiasingGPILSelection;
  return fWrappedProcessAlongStepGPIL;

}

G4VParticleChange* G4BiasingProcessInterface::AlongStepDoIt(const G4Track& track,
							    const G4Step&   step)
{
  // ---------------------------------------
  // -- case outside of volume with biasing:
  // ---------------------------------------
  if ( fCurrentBiasingOperator == 0 )
    {
      if ( fWrappedProcessIsAlong ) return fWrappedProcess->AlongStepDoIt(track, step);
      else
	{
	  fDummyParticleChange->Initialize( track );
	  return fDummyParticleChange;
	}
    }
  
  // -----------------------------------
  // -- case inside volume with biasing:
  // -----------------------------------
  if ( fWrappedProcessIsAlong ) fOccurenceBiasingParticleChange->SetWrappedParticleChange( fWrappedProcess->AlongStepDoIt(track, step) );
  else  
    {
      fOccurenceBiasingParticleChange->SetWrappedParticleChange ( 0 );
      fOccurenceBiasingParticleChange->ProposeTrackStatus( track.GetTrackStatus() );
    }
  G4double weightForNonInteraction (1.0);
  if ( fBiasingInteractionLaw != 0 ) 
    {
      weightForNonInteraction =
  	fPhysicalInteractionLaw->ComputeNonInteractionProbabilityAt(step.GetStepLength()) /
  	fBiasingInteractionLaw ->ComputeNonInteractionProbabilityAt(step.GetStepLength());
      
      fOccurenceBiasingOperation->AlongMoveBy( this, &step, weightForNonInteraction );

      if ( weightForNonInteraction <= 0. )
	{
	  G4ExceptionDescription ed;
	  ed << " Negative non interaction weight : w_NI = " << weightForNonInteraction <<
	    " p_NI(phys) = " <<  fPhysicalInteractionLaw->ComputeNonInteractionProbabilityAt(step.GetStepLength()) <<
	    " p_NI(bias) = " <<  fBiasingInteractionLaw ->ComputeNonInteractionProbabilityAt(step.GetStepLength()) <<
	    " step length = "  <<  step.GetStepLength() <<
	    " biasing interaction law = `" << fBiasingInteractionLaw->GetName() << "'" << G4endl;
	  G4Exception(" G4BiasingProcessInterface::AlongStepDoIt(...)",
		      "BIAS.GEN.04",
		      JustWarning,
		      ed);
	}
      
    }
  
  fOccurenceBiasingParticleChange->SetOccurenceWeightForNonInteraction( weightForNonInteraction );
  
  return fOccurenceBiasingParticleChange;

}

// -- AtRest methods
G4double           G4BiasingProcessInterface::AtRestGetPhysicalInteractionLength(const G4Track&    track,
										 G4ForceCondition* condition)
{
  return  fWrappedProcess->AtRestGetPhysicalInteractionLength(track, condition);
}
G4VParticleChange* G4BiasingProcessInterface::AtRestDoIt(const G4Track& track,
							 const G4Step& step)
{
  return  fWrappedProcess->AtRestDoIt(track, step);
}


G4bool         G4BiasingProcessInterface::IsApplicable(const G4ParticleDefinition& pd)
{
  if ( fWrappedProcess != 0 ) return fWrappedProcess->IsApplicable(pd);
  else                        return true;
}


void       G4BiasingProcessInterface::SetMasterProcess(G4VProcess* masterP)
{
  // -- Master for this process:
  G4VProcess::SetMasterProcess(masterP);
  // -- Master for wrapped process:
  if ( fWrappedProcess != 0 )
    {
      const G4BiasingProcessInterface* thisWrapperMaster = (const G4BiasingProcessInterface *)GetMasterProcess();
      // -- paranoia check:
      G4VProcess* wrappedMaster = 0;
      wrappedMaster = thisWrapperMaster->GetWrappedProcess();
      fWrappedProcess->SetMasterProcess( wrappedMaster );
    }
}
void      G4BiasingProcessInterface::BuildPhysicsTable(const G4ParticleDefinition& pd)
{
  // -- Inform existing operators about start of the run.
  // -- IMPORTANT : as PreparePhysicsTable(...) has been called first for all processes,
  // -- the first/last flags and G4BiasingProcessInterface vector of processes have
  // -- been properly setup.
  if ( fIamFirstGPIL )
    {
      for ( size_t optr = 0 ; optr < ( G4VBiasingOperator::GetBiasingOperators() ).size() ; optr ++)
   	( G4VBiasingOperator::GetBiasingOperators() )[optr]->StartRun( );
    }
  
  if ( fWrappedProcess != 0 )
    {
       fWrappedProcess->BuildPhysicsTable(pd);
    }
}
void    G4BiasingProcessInterface::PreparePhysicsTable(const G4ParticleDefinition& pd)
{
  // -- Let process finding its first/last position in the process manager:
  SetUpFirstLastFlags();
  if ( fWrappedProcess != 0 )
    {
       fWrappedProcess->PreparePhysicsTable(pd);
    }
}
G4bool    G4BiasingProcessInterface::StorePhysicsTable(const G4ParticleDefinition* pd, const G4String& s, G4bool f)
{
  if ( fWrappedProcess != 0 ) return fWrappedProcess->StorePhysicsTable(pd, s, f);
  else                        return false;
}
G4bool G4BiasingProcessInterface::RetrievePhysicsTable(const G4ParticleDefinition* pd, const G4String& s, G4bool f)
{
  if ( fWrappedProcess != 0 ) return fWrappedProcess->RetrievePhysicsTable(pd, s, f);
  else                        return false;
}

void G4BiasingProcessInterface::SetProcessManager(const G4ProcessManager* mgr)
{
  if ( fWrappedProcess != 0 ) fWrappedProcess->SetProcessManager(mgr);
  else                        G4VProcess::SetProcessManager(mgr);
  (fManagerInterfaceMap[mgr]).push_back(this);
  fCoInterfaces = &(fManagerInterfaceMap[mgr]);
  fProcessManager = mgr;
}

const G4ProcessManager* G4BiasingProcessInterface::GetProcessManager()
{
  if ( fWrappedProcess != 0 ) return fWrappedProcess->GetProcessManager();
  else                        return G4VProcess::GetProcessManager();
}
void G4BiasingProcessInterface::BuildWorkerPhysicsTable(const G4ParticleDefinition& pd)
{
  // -- Inform existing operators about start of the run.
  // -- IMPORTANT : as PreparePhysicsTable(...) has been called first for all processes,
  // -- the first/last flags and G4BiasingProcessInterface vector of processes have
  // -- been properly setup.
  if ( fIamFirstGPIL )
    {
      for ( size_t optr = 0 ; optr < ( G4VBiasingOperator::GetBiasingOperators() ).size() ; optr ++)
   	( G4VBiasingOperator::GetBiasingOperators() )[optr]->StartRun( );
    }

  if ( fWrappedProcess != 0 )
    {
      fWrappedProcess->BuildWorkerPhysicsTable(pd);
    }
}
void G4BiasingProcessInterface::PrepareWorkerPhysicsTable(const G4ParticleDefinition& pd)
{
 // -- Let process finding its first/last position in the process manager:
  SetUpFirstLastFlags();
  if ( fWrappedProcess != 0 )
    {
      fWrappedProcess->PrepareWorkerPhysicsTable(pd);
    }  
}
void G4BiasingProcessInterface::ResetNumberOfInteractionLengthLeft()
{
  if ( fWrappedProcess != 0 ) fWrappedProcess->ResetNumberOfInteractionLengthLeft();
}


G4bool G4BiasingProcessInterface::GetIsFirstPostStepGPILInterface( G4bool physOnly ) const
{
  G4int iPhys = ( physOnly ) ? 1 : 0;
  return fFirstLastFlags[IdxFirstLast( 1, 1, iPhys)];
}
G4bool  G4BiasingProcessInterface::GetIsLastPostStepGPILInterface( G4bool physOnly ) const
{
  G4int iPhys = ( physOnly ) ? 1 : 0;
  return fFirstLastFlags[IdxFirstLast( 0, 1, iPhys)];
}
G4bool G4BiasingProcessInterface::GetIsFirstPostStepDoItInterface( G4bool physOnly ) const
{
  G4int iPhys = ( physOnly ) ? 1 : 0;
  return fFirstLastFlags[IdxFirstLast( 1, 0, iPhys)];
}
G4bool  G4BiasingProcessInterface::GetIsLastPostStepDoItInterface( G4bool physOnly ) const
{
  G4int iPhys = ( physOnly ) ? 1 : 0;
  return fFirstLastFlags[IdxFirstLast( 0, 0, iPhys)];
}


G4bool G4BiasingProcessInterface::IsFirstPostStepGPILInterface(G4bool physOnly) const
{
  G4bool isFirst = true;
  const G4ProcessVector* pv = fProcessManager->GetPostStepProcessVector(typeGPIL);
  G4int thisIdx(-1);
  for (G4int i = 0; i < pv->size(); i++ ) if ( (*pv)(i) == this ) { thisIdx = i; break; }
  for ( size_t i = 0; i < fCoInterfaces->size(); i++ )
    {
      if ( ( (*fCoInterfaces)[i]->GetWrappedProcess() != 0 ) || !physOnly )
	{
	  G4int thatIdx(-1);
	  for (G4int j = 0; j < pv->size(); j++ ) if ( (*pv)(j) ==  (*fCoInterfaces)[i] ) { thatIdx = j; break; }
	  if ( thisIdx >  thatIdx )
	    {
	      isFirst = false;
	      break;
	    }
	}
    }
  return isFirst;
}

G4bool G4BiasingProcessInterface::IsLastPostStepGPILInterface(G4bool physOnly) const
{
  G4bool isLast = true;
  const G4ProcessVector* pv = fProcessManager->GetPostStepProcessVector(typeGPIL);
  G4int thisIdx(-1);
  for (G4int i = 0; i < pv->size(); i++ ) if ( (*pv)(i) == this ) { thisIdx = i; break; }
  for ( size_t i = 0; i < fCoInterfaces->size(); i++ )
    {
      if ( ( (*fCoInterfaces)[i]->GetWrappedProcess() != 0  ) || !physOnly )
	{
	  G4int thatIdx(-1);
	  for (G4int j = 0; j < pv->size(); j++ ) if ( (*pv)(j) ==  (*fCoInterfaces)[i] ) { thatIdx = j; break; }
	  if ( thisIdx <  thatIdx )
	    {
	      isLast = false;
	      break;
	    }
	}
    }
  return isLast;  
}

G4bool G4BiasingProcessInterface::IsFirstPostStepDoItInterface(G4bool physOnly) const
{
  G4bool isFirst = true;
  const G4ProcessVector* pv = fProcessManager->GetPostStepProcessVector(typeDoIt);
  G4int thisIdx(-1);
  for (G4int i = 0; i < pv->size(); i++ ) if ( (*pv)(i) == this ) { thisIdx = i; break; }
  for ( size_t i = 0; i < fCoInterfaces->size(); i++ )
    {
      if ( ( (*fCoInterfaces)[i]->GetWrappedProcess() != 0  ) || !physOnly )
	{
	  G4int thatIdx(-1);
	  for (G4int j = 0; j < pv->size(); j++ ) if ( (*pv)(j) ==  (*fCoInterfaces)[i] ) { thatIdx = j; break; }
	  if ( thisIdx >  thatIdx )
	    {
	      isFirst = false;
	      break;
	    }
	}
    }
  return isFirst;
}

G4bool G4BiasingProcessInterface::IsLastPostStepDoItInterface(G4bool physOnly) const
{
  G4bool isLast = true;
  const G4ProcessVector* pv = fProcessManager->GetPostStepProcessVector(typeDoIt);
  G4int thisIdx(-1);
  for (G4int i = 0; i < pv->size(); i++ ) if ( (*pv)(i) == this ) { thisIdx = i; break; }
  for ( size_t i = 0; i < fCoInterfaces->size(); i++ )
    {
      if ( ( (*fCoInterfaces)[i]->GetWrappedProcess() != 0  ) || !physOnly )
	{
	  G4int thatIdx(-1);
	  for (G4int j = 0; j < pv->size(); j++ ) if ( (*pv)(j) ==  (*fCoInterfaces)[i] ) { thatIdx = j; break; }
	  if ( thisIdx <  thatIdx )
	    {
	      isLast = false;
	      break;
	    }
	}
    }
  return isLast;  
}


void G4BiasingProcessInterface::SetUpFirstLastFlags()
{
  for ( G4int iPhys = 0; iPhys < 2; iPhys++ )
    {
      G4bool physOnly = ( iPhys == 1 );
      fFirstLastFlags[IdxFirstLast( 1, 1, iPhys)] = IsFirstPostStepGPILInterface(physOnly);
      fFirstLastFlags[IdxFirstLast( 0, 1, iPhys)] =  IsLastPostStepGPILInterface(physOnly);
      fFirstLastFlags[IdxFirstLast( 1, 0, iPhys)] = IsFirstPostStepDoItInterface(physOnly);
      fFirstLastFlags[IdxFirstLast( 0, 0, iPhys)] =  IsLastPostStepDoItInterface(physOnly);
    }
  
  // -- for itself, for optimization:
  fIamFirstGPIL =  GetIsFirstPostStepGPILInterface( false );
}



void G4BiasingProcessInterface::ResetForUnbiasedTracking()
{
  fOccurenceBiasingOperation  = 0;
  fFinalStateBiasingOperation = 0;
  fNonPhysicsBiasingOperation = 0;
  fBiasingInteractionLaw      = 0;
}
