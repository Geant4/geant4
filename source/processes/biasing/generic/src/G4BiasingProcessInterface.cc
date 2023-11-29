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
#include "G4ParticleChangeForNothing.hh"
#include "G4VBiasingInteractionLaw.hh"
#include "G4InteractionLawPhysical.hh"
#include "G4ProcessManager.hh"
#include "G4BiasingAppliedCase.hh"
#include "G4ParallelGeometriesLimiterProcess.hh"

G4Cache<G4bool>                     G4BiasingProcessInterface::fResetInteractionLaws;// = true;
G4Cache<G4bool>                     G4BiasingProcessInterface::fCommonStart;//          = true;
G4Cache<G4bool>                     G4BiasingProcessInterface::fCommonEnd;//            = true;
G4Cache<G4bool>                     G4BiasingProcessInterface::fDoCommonConfigure;

G4BiasingProcessInterface::G4BiasingProcessInterface(G4String name)
  :  G4VProcess                           ( name    ),
     fCurrentTrack                        ( nullptr ),
     fPreviousStepSize (-1.0), fCurrentMinimumStep( -1.0 ), fProposedSafety ( -1.0),
             fOccurenceBiasingOperation( nullptr ),         fFinalStateBiasingOperation( nullptr ),         fNonPhysicsBiasingOperation( nullptr ),
     fPreviousOccurenceBiasingOperation( nullptr ), fPreviousFinalStateBiasingOperation( nullptr ), fPreviousNonPhysicsBiasingOperation( nullptr ),
     fResetWrappedProcessInteractionLength( true    ),
     fWrappedProcess                      ( nullptr ),
     fIsPhysicsBasedBiasing               ( false   ),
     fWrappedProcessIsAtRest              ( false   ),
     fWrappedProcessIsAlong               ( false   ),
     fWrappedProcessIsPost                ( false   ),
     fWrappedProcessPostStepGPIL          ( -1.0    ),
     fBiasingPostStepGPIL                 ( -1.0    ),
     fWrappedProcessInteractionLength     ( -1.0    ),
     fWrappedProcessForceCondition        ( NotForced ),
     fBiasingForceCondition               ( NotForced ),
     fWrappedProcessAlongStepGPIL         ( -1.0    ),
     fBiasingAlongStepGPIL                ( -1.0    ),
     fWrappedProcessGPILSelection         ( NotCandidateForSelection ),
     fBiasingGPILSelection                ( NotCandidateForSelection ),
     fBiasingInteractionLaw               ( nullptr ),
     fPreviousBiasingInteractionLaw       ( nullptr ),
     fPhysicalInteractionLaw              ( nullptr ),
     fOccurenceBiasingParticleChange      ( nullptr ),
     fDummyParticleChange                 ( nullptr ),
     fIamFirstGPIL                        ( false   ),
     fProcessManager                      ( nullptr ),
     fSharedData                          ( nullptr )
{
  for (G4int i = 0 ; i < 8 ; i++)  fFirstLastFlags[i] = false;
  fResetInteractionLaws.Put( true );
  fCommonStart         .Put( true );
  fCommonEnd           .Put( true );
  fDoCommonConfigure   .Put( true );
}


G4BiasingProcessInterface::G4BiasingProcessInterface(G4VProcess* wrappedProcess,
						     G4bool wrappedIsAtRest, G4bool wrappedIsAlongStep, G4bool wrappedIsPostStep,
						     G4String useThisName)
  : G4VProcess( useThisName != "" ? useThisName : "biasWrapper("+wrappedProcess->GetProcessName()+")",
	       wrappedProcess->GetProcessType()),
    fCurrentTrack                         ( nullptr            ),
    fPreviousStepSize (-1.0), fCurrentMinimumStep( -1.0 ), fProposedSafety ( -1.0),
            fOccurenceBiasingOperation( nullptr ),         fFinalStateBiasingOperation( nullptr ),         fNonPhysicsBiasingOperation( nullptr ),
    fPreviousOccurenceBiasingOperation( nullptr ), fPreviousFinalStateBiasingOperation( nullptr ), fPreviousNonPhysicsBiasingOperation( nullptr ),
    fResetWrappedProcessInteractionLength ( false              ),
    fWrappedProcess                       ( wrappedProcess     ),
    fIsPhysicsBasedBiasing                ( true               ),
    fWrappedProcessIsAtRest               ( wrappedIsAtRest    ),
    fWrappedProcessIsAlong                ( wrappedIsAlongStep ),
    fWrappedProcessIsPost                 ( wrappedIsPostStep  ),
    fWrappedProcessPostStepGPIL           ( -1.0               ),
    fBiasingPostStepGPIL                  ( -1.0               ),
    fWrappedProcessInteractionLength      ( -1.0               ),
    fWrappedProcessForceCondition         ( NotForced          ),
    fBiasingForceCondition                ( NotForced          ),
    fWrappedProcessAlongStepGPIL          ( -1.0               ),
    fBiasingAlongStepGPIL                 ( -1.0               ),
    fWrappedProcessGPILSelection          ( NotCandidateForSelection ),
    fBiasingGPILSelection                 ( NotCandidateForSelection ),
    fBiasingInteractionLaw                ( nullptr            ),
    fPreviousBiasingInteractionLaw        ( nullptr            ),
    fPhysicalInteractionLaw               ( nullptr            ),
    fOccurenceBiasingParticleChange       ( nullptr            ),
    fIamFirstGPIL                         ( false              ),
    fProcessManager                       ( nullptr            ),
    fSharedData                           ( nullptr            )
{
  for (G4int i = 0 ; i < 8 ; i++)  fFirstLastFlags[i] = false;
  fResetInteractionLaws.Put( true );
  fCommonStart.Put(true);
  fCommonEnd.Put(true);
  fDoCommonConfigure.Put(true);
  
  SetProcessSubType(fWrappedProcess->GetProcessSubType());

  // -- create physical interaction law:
  fPhysicalInteractionLaw         = new            G4InteractionLawPhysical("PhysicalInteractionLawFor("+GetProcessName()+")");
  // -- instantiate particle change wrapper for occurrence biaising:
  fOccurenceBiasingParticleChange = new G4ParticleChangeForOccurenceBiasing("biasingPCfor"+GetProcessName());
  // -- instantiate a "do nothing" particle change:
  fDummyParticleChange            = new          G4ParticleChangeForNothing();
}



G4BiasingProcessInterface::~G4BiasingProcessInterface()
{
  if ( fPhysicalInteractionLaw  != 0   ) delete fPhysicalInteractionLaw;
  if ( fOccurenceBiasingParticleChange ) delete fOccurenceBiasingParticleChange;
  if ( fDummyParticleChange            ) delete fDummyParticleChange;
}


const G4BiasingProcessSharedData* G4BiasingProcessInterface::GetSharedData( const G4ProcessManager* mgr )
{
  G4MapCache< const G4ProcessManager*, 
	      G4BiasingProcessSharedData* >::const_iterator itr =   G4BiasingProcessSharedData::fSharedDataMap.Find( mgr );
  if ( itr !=  G4BiasingProcessSharedData::fSharedDataMap.End( ) )
    {
      return (*itr).second;
    }
  else return 0;
}


void G4BiasingProcessInterface::StartTracking(G4Track* track)
{
  fCurrentTrack = track;
  if ( fIsPhysicsBasedBiasing ) fWrappedProcess->StartTracking(fCurrentTrack);
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
      fCommonEnd  .Put( true  );// = true;
      
      fSharedData-> fCurrentBiasingOperator = 0;
      fSharedData->fPreviousBiasingOperator = 0;

      // -- §§ Add a  "fSharedData->nStarting" here and outside bracket "fSharedData->nStarting++" and " if (fSharedData->nStarting) == fSharedData->(vector interface length)"
      // -- §§ call to the loop "StartTracking" of operators"
      
      for ( size_t optr = 0 ; optr < ( G4VBiasingOperator::GetBiasingOperators() ).size() ; optr ++)
	( G4VBiasingOperator::GetBiasingOperators() )[optr]->StartTracking( fCurrentTrack );
    }
}


void G4BiasingProcessInterface::EndTracking()
{
  if ( fIsPhysicsBasedBiasing ) fWrappedProcess->EndTracking();
  if ( fSharedData->fCurrentBiasingOperator) (fSharedData->fCurrentBiasingOperator)->ExitingBiasing( fCurrentTrack, this ); 
  fBiasingInteractionLaw = 0;

  // -- Inform operators of end of tracking:
  if ( fCommonEnd.Get() )
    {
      fCommonEnd  .Put( false );// = false;
      fCommonStart.Put( true  );// = true;
      
      for ( size_t optr = 0 ; optr < ( G4VBiasingOperator::GetBiasingOperators() ).size() ; optr ++)
	( G4VBiasingOperator::GetBiasingOperators() )[optr]->EndTracking( );

      // -- §§ for above loop, do as in StartTracking.
    }
}



G4double G4BiasingProcessInterface::PostStepGetPhysicalInteractionLength( const G4Track&                track,
									  G4double           previousStepSize,
									  G4ForceCondition*         condition )
{

  // ---------------------------------------------------------------------------------------------------
  // -- The "biasing process master" takes care of updating the biasing operator, and for all biasing
  // -- processes it invokes the PostStepGPIL of physical wrapped processes (anticipate stepping manager
  // -- call ! ) to make all cross-sections updated with current step, and hence available before the
  // -- first call to the biasing operator.
  // ---------------------------------------------------------------------------------------------------
  if ( fIamFirstGPIL )
    {
      // -- Update previous biasing operator, and assume the operator stays the same by
      // -- default and that it is not left at the beginning of this step. These
      // -- assumptions might be wrong if there is a volume change (in paralllel or
      // -- mass geometries) in what case the flags will be updated.
      fSharedData->fPreviousBiasingOperator = fSharedData->fCurrentBiasingOperator;
      fSharedData->fIsNewOperator           = false;
      fSharedData->fLeavingPreviousOperator = false;
      // -- If new volume, either in mass or parallel geometries, get possible new biasing operator:
      // -------------------------------------------------------------------------------------------
      // -- Get biasing operator in parallel geometries:
      G4bool firstStepInParallelVolume = false;
      if ( fSharedData->fParallelGeometriesLimiterProcess )
	{
	  G4VBiasingOperator* newParallelOperator( nullptr );
	  G4bool firstStep = ( track.GetCurrentStepNumber() == 1 );
	  size_t iParallel = 0;
	  for ( auto wasLimiting : fSharedData->fParallelGeometriesLimiterProcess->GetWasLimiting() )
	    {
	      if ( firstStep || wasLimiting )
		{
		  firstStepInParallelVolume = true;
		  
		  auto tmpParallelOperator = G4VBiasingOperator::GetBiasingOperator( (fSharedData->fParallelGeometriesLimiterProcess->GetCurrentVolumes()[iParallel])
										     ->GetLogicalVolume()                                                             );
		  if ( newParallelOperator )
		    {
		      if ( tmpParallelOperator )
			{
			  G4ExceptionDescription ed;
			  ed << " Several biasing operators are defined at the same place in parallel geometries ! Found:\n";
			  ed << "    - `" << newParallelOperator->GetName() << "' and \n";
			  ed << "    - `" << tmpParallelOperator->GetName() << "'.\n";
			  ed << " Keeping `" << newParallelOperator->GetName() << "'. Behavior not guaranteed ! Please consider having only one operator at a place. " << G4endl;
			  G4Exception(" G4BiasingProcessInterface::PostStepGetPhysicalInteractionLength(...)",
				      "BIAS.GEN.30",
				      JustWarning,
				      ed);
			}
		    }
		  else newParallelOperator = tmpParallelOperator;
		}
	      iParallel++;
	    }
	  fSharedData->fParallelGeometryOperator = newParallelOperator;
	} // -- end of " if ( fSharedData->fParallelGeometriesLimiterProcess )"

      // -- Get biasing operator in mass geometry:
      // -- [§§ Note : bug with this first step ? Does not work if previous step was concurrently limited with geometry. Might make use of safety at last point ?]
      G4bool  firstStepInVolume = ( (track.GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) || (track.GetCurrentStepNumber() == 1) );
      //      fSharedData->fIsNewOperator           = false;
      //      fSharedData->fLeavingPreviousOperator = false;
      if ( firstStepInVolume )
	{
	  G4VBiasingOperator* newOperator = G4VBiasingOperator::GetBiasingOperator( track.GetVolume()->GetLogicalVolume() );
	  fSharedData->fMassGeometryOperator = newOperator;
	  if ( ( newOperator != nullptr ) && ( fSharedData->fParallelGeometryOperator != nullptr ) )
	    {
	      G4ExceptionDescription ed;
	      ed << " Biasing operators are defined at the same place in mass and parallel geometries ! Found:\n";
	      ed << "    - `" << fSharedData->fParallelGeometryOperator->GetName() << "' in parallel geometry and \n";
	      ed << "    - `" << newOperator->GetName() << "' in mass geometry.\n";
	      ed << " Keeping `" << fSharedData->fParallelGeometryOperator->GetName() << "'. Behavior not guaranteed ! Please consider having only one operator at a place. " << G4endl;
	      G4Exception(" G4BiasingProcessInterface::PostStepGetPhysicalInteractionLength(...)",
			  "BIAS.GEN.31",
			  JustWarning,
			  ed);
	    }
	}

      // -- conclude the operator selection, giving priority to parallel geometry (as told in exception message BIAS.GEN.30):
      if ( firstStepInVolume || firstStepInParallelVolume )
	{
	  G4VBiasingOperator*           newOperator = fSharedData->fParallelGeometryOperator;
	  if ( newOperator == nullptr ) newOperator = fSharedData->fMassGeometryOperator;
	  
	  fSharedData->fCurrentBiasingOperator = newOperator ;

	  if ( newOperator != fSharedData->fPreviousBiasingOperator )
	    {
	      fSharedData->fLeavingPreviousOperator = ( fSharedData->fPreviousBiasingOperator != nullptr ) ;
	      fSharedData->fIsNewOperator           = ( newOperator != nullptr );
	    }
	}
      

      // -- calls to wrapped process PostStepGPIL's:
      // -------------------------------------------
      // -- Each physics wrapper process has its
      // --   fWrappedProcessPostStepGPIL      ,
      // --   fWrappedProcessForceCondition    ,
      // --   fWrappedProcessInteractionLength
      // -- updated.
      if ( fSharedData->fCurrentBiasingOperator != nullptr )
	{
	  for ( size_t i = 0 ; i < (fSharedData->fPhysicsBiasingProcessInterfaces).size(); i++ )
	    (fSharedData->fPhysicsBiasingProcessInterfaces)[i]->InvokeWrappedProcessPostStepGPIL( track, previousStepSize, condition );
	}
    } // -- end of "if ( fIamFirstGPIL )"

  

  // -- Remember previous operator and proposed operations, if any, and reset:
  // -------------------------------------------------------------------------
  // -- remember only in case some biasing might be called
  if ( ( fSharedData->fPreviousBiasingOperator != 0 ) ||
       ( fSharedData->fCurrentBiasingOperator  != 0 )    )
    {
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
      // fWrappedProcessPostStepGPIL      : updated by InvokeWrappedProcessPostStepGPIL(...) above
      fBiasingPostStepGPIL                =                     DBL_MAX;
      // fWrappedProcessInteractionLength : updated by InvokeWrappedProcessPostStepGPIL(...) above; inverse of analog cross-section.
      // fWrappedProcessForceCondition    : updated by InvokeWrappedProcessPostStepGPIL(...) above
      fBiasingForceCondition              =                   NotForced;
      fWrappedProcessAlongStepGPIL        =                     DBL_MAX;
      fBiasingAlongStepGPIL               =                     DBL_MAX;
      fWrappedProcessGPILSelection        =    NotCandidateForSelection;
      fBiasingGPILSelection               =    NotCandidateForSelection;
      // -- for helper:
      fPreviousStepSize                   =            previousStepSize;
    }

  
  // -- previous step size value; it is switched to zero if resetting a wrapped process:
  // -- (same trick used than in InvokedWrappedProcessPostStepGPIL )
  G4double usedPreviousStepSize = previousStepSize;

  // ----------------------------------------------
  // -- If leaving a biasing operator, let it know:
  // ----------------------------------------------
  if ( fSharedData->fLeavingPreviousOperator )
    {
      (fSharedData->fPreviousBiasingOperator)->ExitingBiasing( &track, this ); 
      // -- if no further biasing operator, reset process behavior to standard tracking:
      if ( fSharedData->fCurrentBiasingOperator == 0 )
	{
	  ResetForUnbiasedTracking();
	  if ( fIsPhysicsBasedBiasing )
	    {
	      // -- if the physics process has been under occurrence biasing, reset it:
	      if ( fResetWrappedProcessInteractionLength )
		{
		  fResetWrappedProcessInteractionLength = false;
		  fWrappedProcess->ResetNumberOfInteractionLengthLeft();
		  // -- We set "previous step size" as 0.0, to let the process believe this is first step:
		  usedPreviousStepSize = 0.0;
		}
	    }
	}
    }
  

  // --------------------------------------------------------------
  // -- no operator : analog tracking if physics-based, or nothing:
  // --------------------------------------------------------------
  if ( fSharedData->fCurrentBiasingOperator == 0 )
    {
      // -- take note of the "usedPreviousStepSize" value:
      if ( fIsPhysicsBasedBiasing ) return fWrappedProcess->PostStepGetPhysicalInteractionLength(track, usedPreviousStepSize, condition);
      else
  	{
  	  *condition = NotForced;
  	  return DBL_MAX;
  	}
    }

  // --------------------------------------------------
  // -- A biasing operator exists. Proceed with
  // -- treating non-physics and physics biasing cases:
  //---------------------------------------------------
  
  // -- non-physics-based biasing case:
  // ----------------------------------
  if ( !fIsPhysicsBasedBiasing )
    {  
      fNonPhysicsBiasingOperation = (fSharedData->fCurrentBiasingOperator)->GetProposedNonPhysicsBiasingOperation( &track, this );
      if ( fNonPhysicsBiasingOperation == 0 )
	{
	  *condition = NotForced;
	  return DBL_MAX;
	}
      return fNonPhysicsBiasingOperation->DistanceToApplyOperation(&track, previousStepSize, condition);
    }
  

  // -- Physics based biasing case:
  // ------------------------------
  // -- Ask for possible GPIL biasing operation:
  fOccurenceBiasingOperation = (fSharedData->fCurrentBiasingOperator)->GetProposedOccurenceBiasingOperation( &track, this );


  // -- no operation for occurrence biasing, analog GPIL returns the wrapped process GPIL and condition values
  if ( fOccurenceBiasingOperation == 0 )
    {
      *condition = fWrappedProcessForceCondition;
      return       fWrappedProcessPostStepGPIL;
    }

  // -- A valid GPIL biasing operation has been proposed:
  // -- 0) remember wrapped process will need to be reset on biasing exit, if particle survives:
  fResetWrappedProcessInteractionLength = true;
  // -- 1) update process interaction length for reference analog interaction law ( fWrappedProcessInteractionLength updated/collected above):
  fPhysicalInteractionLaw->SetPhysicalCrossSection( 1.0 / fWrappedProcessInteractionLength );
  // -- 2) Collect biasing interaction law:
  // --    The interaction law pointer is collected as a const pointer to the interaction law object.
  // --    This interaction law will be kept under control of the biasing operation, which is the only
  // --    entity that will change the state of the biasing interaction law.
  // --    The force condition for biasing is asked at the same time, passing the analog one as default:
  fBiasingForceCondition           = fWrappedProcessForceCondition;
  fBiasingInteractionLaw           = fOccurenceBiasingOperation->ProvideOccurenceBiasingInteractionLaw( this, fBiasingForceCondition );
  // -- 3) Ask operation to sample the biasing interaction law:
  fBiasingPostStepGPIL             = fBiasingInteractionLaw->GetSampledInteractionLength();

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
  if ( fSharedData->fCurrentBiasingOperator == 0 ) return fWrappedProcess->PostStepDoIt(track, step);
  
  // ----------------------------
  // -- non-physics biasing case:
  // ----------------------------
  if ( !fIsPhysicsBasedBiasing )
    {
      G4VParticleChange* particleChange = fNonPhysicsBiasingOperation->GenerateBiasingFinalState( &track, &step );
      (fSharedData->fCurrentBiasingOperator)->ReportOperationApplied( this, BAC_NonPhysics, fNonPhysicsBiasingOperation, particleChange );
      return particleChange;
    }

  // -- physics biasing case:
  // ------------------------
  // -- It proceeds with the following logic:
  // -- 1) Obtain the final state
  // --    This final state may be analog or biased.
  // --    The biased final state is obtained through a biasing operator
  // --    returned by the operator.
  // -- 2) The biased final state may be asked to be "force as it is"
  // --    in what case the particle change is returned as is to the
  // --    stepping.
  // --    In all other cases (analog final state or biased final but
  // --    not forced) the final state weight may be modified by the
  // --    occurrence biasing, if such an occurrence biasing is at play.
  
  // -- Get final state, biased or analog:
  G4VParticleChange*   finalStateParticleChange;
  G4BiasingAppliedCase BAC;
  fFinalStateBiasingOperation = (fSharedData->fCurrentBiasingOperator)->GetProposedFinalStateBiasingOperation( &track, this );
  // -- Flag below is to force the biased generated particle change to be returned "as is" to the stepping, disregarding there
  // -- was or not a occurrence biasing that would apply. Weight relevance under full responsibility of the biasing operation.
  G4bool forceBiasedFinalState = false;
  if ( fFinalStateBiasingOperation != 0 )
    {
      finalStateParticleChange = fFinalStateBiasingOperation->ApplyFinalStateBiasing( this, &track, &step, forceBiasedFinalState );
      BAC = BAC_FinalState;
    }
  else
    {
      finalStateParticleChange = fWrappedProcess->PostStepDoIt(track, step);
      BAC =  BAC_None ;
    }
  
  // -- if no occurrence biasing operation, we're done:
  if ( fOccurenceBiasingOperation == 0 )
    {
      (fSharedData->fCurrentBiasingOperator)->ReportOperationApplied( this, BAC, fFinalStateBiasingOperation, finalStateParticleChange );
      return finalStateParticleChange;
    }
  
  // -- if biased final state has been asked to be forced, we're done:
  if ( forceBiasedFinalState )
    {
      (fSharedData->fCurrentBiasingOperator)->ReportOperationApplied( this, BAC, fFinalStateBiasingOperation, finalStateParticleChange );
      return finalStateParticleChange;
    }
  

  // -- If occurrence biasing, applies the occurrence biasing weight correction on top of final state (biased or not):
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
	" XS_I(phys) = "       << fBiasingInteractionLaw ->ComputeEffectiveCrossSectionAt(step.GetStepLength()) <<
	" XS_I(bias) = "       << fPhysicalInteractionLaw->ComputeEffectiveCrossSectionAt(step.GetStepLength()) <<
	" step length = "      << step.GetStepLength() <<
	" Interaction law = `" << fBiasingInteractionLaw << "'" <<
	G4endl;
      G4Exception(" G4BiasingProcessInterface::PostStepDoIt(...)",
		  "BIAS.GEN.03",
		  JustWarning,
		  ed);
      
    }
  
  (fSharedData->fCurrentBiasingOperator)->ReportOperationApplied( this,                        BAC,
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
  if ( fSharedData->fCurrentBiasingOperator == 0 )
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

  
  // -----------------------------------------------------------
  // -- From here we have an valid occurrence biasing operation:
  // -----------------------------------------------------------
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
  if ( fSharedData->fCurrentBiasingOperator == 0 )
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
      // -- paranoia check: (?)
      G4VProcess* wrappedMaster = 0;
      wrappedMaster = thisWrapperMaster->GetWrappedProcess();
      fWrappedProcess->SetMasterProcess( wrappedMaster );
    }
}


void      G4BiasingProcessInterface::BuildPhysicsTable(const G4ParticleDefinition& pd)
{
  // -- Sequential mode : called second (after PreparePhysicsTable(..))
  // -- MT mode         : called second (after PreparePhysicsTable(..)) by master thread.
  // --                   Corresponding process instance not used then by tracking.
  // -- PreparePhysicsTable(...) has been called first for all processes,
  // -- so the first/last flags and G4BiasingProcessInterface vector of processes have
  // -- been properly setup, fIamFirstGPIL is valid.
  if ( fWrappedProcess != 0 )
    {
      fWrappedProcess->BuildPhysicsTable(pd);
    }

  if ( fIamFirstGPIL )
    {
      // -- Re-order vector of processes to match that of the GPIL
      // -- (made for fIamFirstGPIL, but important is to have it made once):
      ReorderBiasingVectorAsGPIL();
      // -- Let operators to configure themselves for the master thread or for sequential mode.
      // -- Intended here is in particular the registration to physics model catalog.
      // -- The fDoCommonConfigure is to ensure that this Configure is made by only one process (othewise each first process makes the call):
      if ( fDoCommonConfigure.Get() )
	{
	  for ( size_t optr = 0 ; optr < ( G4VBiasingOperator::GetBiasingOperators() ).size() ; optr ++)
	    ( G4VBiasingOperator::GetBiasingOperators() )[optr]->Configure( );
	  fDoCommonConfigure.Put(false);
	}
      
    }
}


void    G4BiasingProcessInterface::PreparePhysicsTable(const G4ParticleDefinition& pd)
{
  // -- Sequential mode : called first (before BuildPhysicsTable(..))
  // -- MT mode         : called first (before BuildPhysicsTable(..)) by master thread.
  // --                   Corresponding process instance not used then by tracking.
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

  // -- initialize fSharedData pointer:
  if (  G4BiasingProcessSharedData::fSharedDataMap.Find(mgr) == G4BiasingProcessSharedData::fSharedDataMap.End() )
    {
      fSharedData =  new G4BiasingProcessSharedData( mgr );
      G4BiasingProcessSharedData::fSharedDataMap[mgr] = fSharedData;
    }
  else fSharedData =  G4BiasingProcessSharedData::fSharedDataMap[mgr] ;
  // -- augment list of co-operating processes:
  fSharedData->       fBiasingProcessInterfaces.push_back( this );
  fSharedData-> fPublicBiasingProcessInterfaces.push_back( this );
  if ( fIsPhysicsBasedBiasing ) 
    {
      fSharedData->       fPhysicsBiasingProcessInterfaces.push_back( this );
      fSharedData-> fPublicPhysicsBiasingProcessInterfaces.push_back( this );
    }
  else
    {
      fSharedData->       fNonPhysicsBiasingProcessInterfaces.push_back( this );
      fSharedData-> fPublicNonPhysicsBiasingProcessInterfaces.push_back( this );
    }
  // -- remember process manager:
  fProcessManager = mgr;
}


const G4ProcessManager* G4BiasingProcessInterface::GetProcessManager()
{
  if ( fWrappedProcess != 0 ) return fWrappedProcess->GetProcessManager();
  else                        return G4VProcess::GetProcessManager();
}


void G4BiasingProcessInterface::BuildWorkerPhysicsTable(const G4ParticleDefinition& pd)
{
  // -- Sequential mode : not called
  // -- MT mode         : called after PrepareWorkerPhysicsTable(..)
  // -- PrepareWorkerPhysicsTable(...) has been called first for all processes,
  // -- so the first/last flags and G4BiasingProcessInterface vector of processes have
  // -- been properly setup, fIamFirstGPIL is valid.
  if ( fWrappedProcess != 0 )
    {
      fWrappedProcess->BuildWorkerPhysicsTable(pd);
    }

  if ( fIamFirstGPIL )
    {
      // -- Re-order vector of processes to match that of the GPIL
      // -- (made for fIamFirstGPIL, but important is to have it made once):
      ReorderBiasingVectorAsGPIL();
      // -- Let operators to configure themselves for the worker thread, if needed.
      // -- Registration to physics model catalog **IS NOT** to be made here, but in Configure().
      // -- The fDoCommonConfigure is to ensure that this Configure is made by only one process (othewise each first process makes the call):
      if ( fDoCommonConfigure.Get() )
	{
	  for ( size_t optr = 0 ; optr < ( G4VBiasingOperator::GetBiasingOperators() ).size() ; optr ++)
	    ( G4VBiasingOperator::GetBiasingOperators() )[optr]->ConfigureForWorker( );
	  fDoCommonConfigure.Put(false);
	}
    }
}


void G4BiasingProcessInterface::PrepareWorkerPhysicsTable(const G4ParticleDefinition& pd)
{
  // -- Sequential mode : not called
  // -- MT mode         : called first, before BuildWorkerPhysicsTable(..)
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
  for (G4int i = 0; i < (G4int)pv->size(); ++i )
    if ( (*pv)(i) == this ) { thisIdx = i; break; }
  if ( thisIdx < 0 ) return false; // -- to ignore pure along processes
  for ( std::size_t i = 0; i < (fSharedData->fBiasingProcessInterfaces).size(); ++i )
    {
      if ( (fSharedData->fBiasingProcessInterfaces)[i]->fIsPhysicsBasedBiasing || !physOnly )
	{
	  G4int thatIdx(-1);
	  for (G4int j = 0; j < (G4int)pv->size(); ++j )
            if ( (*pv)(j) == (fSharedData->fBiasingProcessInterfaces)[i] )
              { thatIdx = j; break; }
	  if ( thatIdx >= 0 ) // -- to ignore pure along processes
	    {
	      if ( thisIdx >  thatIdx )
		{
		  isFirst = false;
		  break;
		}
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
  for (G4int i = 0; i < (G4int)pv->size(); ++i )
    if ( (*pv)(i) == this ) { thisIdx = i; break; }
  if ( thisIdx < 0 ) return false; // -- to ignore pure along processes
  for ( std::size_t i = 0; i < (fSharedData->fBiasingProcessInterfaces).size(); ++i )
    {
      if ( (fSharedData->fBiasingProcessInterfaces)[i]->fIsPhysicsBasedBiasing || !physOnly )
	{
	  G4int thatIdx(-1);
	  for (G4int j = 0; j < (G4int)pv->size(); ++j )
            if ( (*pv)(j) == (fSharedData->fBiasingProcessInterfaces)[i] )
              { thatIdx = j; break; }
	  if ( thatIdx >= 0 ) // -- to ignore pure along processes
	    {
	      if ( thisIdx <  thatIdx )
		{
		  isLast = false;
		  break;
		}
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
  for (G4int i = 0; i < (G4int)pv->size(); ++i )
    if ( (*pv)(i) == this ) { thisIdx = i; break; }
  if ( thisIdx < 0 ) return false; // -- to ignore pure along processes
  for ( std::size_t i = 0; i < (fSharedData->fBiasingProcessInterfaces).size(); ++i )
    {
      if ( (fSharedData->fBiasingProcessInterfaces)[i]->fIsPhysicsBasedBiasing || !physOnly )
	{
	  G4int thatIdx(-1);
	  for (G4int j = 0; j < (G4int)pv->size(); ++j )
            if ( (*pv)(j) == (fSharedData->fBiasingProcessInterfaces)[i] )
              { thatIdx = j; break; }
	  if ( thatIdx >= 0 ) // -- to ignore pure along processes
	    {
	      if ( thisIdx >  thatIdx )
		{
		  isFirst = false;
		  break;
		}
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
  for (G4int i = 0; i < (G4int)pv->size(); ++i )
    if ( (*pv)(i) == this ) { thisIdx = i; break; }
  if ( thisIdx < 0 ) return false; // -- to ignore pure along processes
  for ( std::size_t i = 0; i < (fSharedData->fBiasingProcessInterfaces).size(); ++i )
    {
      if ( (fSharedData->fBiasingProcessInterfaces)[i]->fIsPhysicsBasedBiasing || !physOnly )
	{
	  G4int thatIdx(-1);
	  for (G4int j = 0; j < (G4int)pv->size(); ++j )
            if ( (*pv)(j) == (fSharedData->fBiasingProcessInterfaces)[i] )
              { thatIdx = j; break; }
	  if ( thatIdx >= 0 ) // -- to ignore pure along processes
	    {
	      if ( thisIdx <  thatIdx )
		{
		  isLast = false;
		  break;
		}
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


void G4BiasingProcessInterface::InvokeWrappedProcessPostStepGPIL( const G4Track&               track,
								  G4double          previousStepSize,
								  G4ForceCondition*        condition )
{
  G4double usedPreviousStepSize = previousStepSize;
  // -- if the physics process has been under occurrence biasing in the previous step
  // -- we reset it, as we don't know if it will be biased again or not in this
  // -- step. The pity is that PostStepGPIL and interaction length (cross-section)
  // -- calculations are done both in the PostStepGPIL of the process, while here we
  // -- are just interested in the calculation of the cross-section. This is a pity
  // -- as this forces to re-generated a random number for nothing.
  if ( fResetWrappedProcessInteractionLength )
    {
      fResetWrappedProcessInteractionLength = false;
      fWrappedProcess->ResetNumberOfInteractionLengthLeft();
      // -- We set "previous step size" as 0.0, to let the process believe this is first step:
      usedPreviousStepSize = 0.0;
    }
  // -- GPIL response:
  fWrappedProcessPostStepGPIL      = fWrappedProcess->PostStepGetPhysicalInteractionLength(track, usedPreviousStepSize, condition);
  fWrappedProcessForceCondition    = *condition;
  // -- and (inverse) cross-section:
  fWrappedProcessInteractionLength = fWrappedProcess->GetCurrentInteractionLength();
}


void G4BiasingProcessInterface::ReorderBiasingVectorAsGPIL()
{
  // -- re-order vector of processes to match that of the GPIL:
  std::vector < G4BiasingProcessInterface* > tmpProcess ( fSharedData->fBiasingProcessInterfaces );
  ( fSharedData -> fBiasingProcessInterfaces                 ) . clear();
  ( fSharedData -> fPhysicsBiasingProcessInterfaces          ) . clear();
  ( fSharedData -> fNonPhysicsBiasingProcessInterfaces       ) . clear();
  ( fSharedData -> fPublicBiasingProcessInterfaces           ) . clear();
  ( fSharedData -> fPublicPhysicsBiasingProcessInterfaces    ) . clear();
  ( fSharedData -> fPublicNonPhysicsBiasingProcessInterfaces ) . clear();
  
  const G4ProcessVector* pv = fProcessManager->GetPostStepProcessVector(typeGPIL);
  for (G4int i = 0; i < (G4int)pv->size(); ++i ) 
    {
      for ( std::size_t j = 0; j < tmpProcess.size(); ++j )
        {
          if ( (*pv)(i) == tmpProcess[j] )
            { 
              ( fSharedData -> fBiasingProcessInterfaces                     ) . push_back( tmpProcess[j] );
              ( fSharedData -> fPublicBiasingProcessInterfaces               ) . push_back( tmpProcess[j] );
              if ( tmpProcess[j] -> fIsPhysicsBasedBiasing )
                {
                  ( fSharedData -> fPhysicsBiasingProcessInterfaces          ) . push_back( tmpProcess[j] );
                  ( fSharedData -> fPublicPhysicsBiasingProcessInterfaces    ) . push_back( tmpProcess[j] ); 
                }
              else
                {
                  ( fSharedData -> fNonPhysicsBiasingProcessInterfaces       ) . push_back( tmpProcess[j] );
                  ( fSharedData -> fPublicNonPhysicsBiasingProcessInterfaces ) . push_back( tmpProcess[j] );  
                }
              break;
            }
        }
    }
}
