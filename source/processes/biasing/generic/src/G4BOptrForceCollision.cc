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
#include "G4BOptrForceCollision.hh"
#include "G4BOptrForceCollisionTrackData.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4PhysicsModelCatalog.hh"

#include "G4BOptnForceCommonTruncatedExp.hh"
#include "G4ILawCommonTruncatedExp.hh"
#include "G4BOptnForceFreeFlight.hh"
#include "G4BOptnCloning.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4SystemOfUnits.hh"

// -- §§ consider calling other constructor, thanks to C++11
G4BOptrForceCollision::G4BOptrForceCollision(G4String particleName, G4String name)
  : G4VBiasingOperator(name),
    fForceCollisionModelID(-1),
    fCurrentTrack(nullptr),
    fCurrentTrackData(nullptr),
    fInitialTrackWeight(-1.0),
    fSetup(true)
{
  fSharedForceInteractionOperation = new G4BOptnForceCommonTruncatedExp("SharedForceInteraction");
  fCloningOperation                = new G4BOptnCloning("Cloning");
  fParticleToBias = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  
  if ( fParticleToBias == 0 )
    {
      G4ExceptionDescription ed;
      ed << " Particle `" << particleName << "' not found !" << G4endl;
      G4Exception(" G4BOptrForceCollision::G4BOptrForceCollision(...)",
		  "BIAS.GEN.07",
		  JustWarning,
		  ed);
    }
}


G4BOptrForceCollision::G4BOptrForceCollision(const G4ParticleDefinition* particle, G4String name)
  : G4VBiasingOperator(name),
    fForceCollisionModelID(-1),
    fCurrentTrack(nullptr),
    fCurrentTrackData(nullptr),
    fInitialTrackWeight(-1.0),
    fSetup(true)
{
  fSharedForceInteractionOperation = new G4BOptnForceCommonTruncatedExp("SharedForceInteraction");
  fCloningOperation                = new G4BOptnCloning("Cloning");
  fParticleToBias                  = particle;
}


G4BOptrForceCollision::~G4BOptrForceCollision()
{
  for ( std::map< const G4BiasingProcessInterface*, G4BOptnForceFreeFlight* >::iterator it = fFreeFlightOperations.begin() ;
	it != fFreeFlightOperations.end() ;
	it++ ) delete (*it).second;
  delete fSharedForceInteractionOperation;
  delete fCloningOperation;
}


void G4BOptrForceCollision::Configure()
{
  // -- Create ID for force collision:
  fForceCollisionModelID = G4PhysicsModelCatalog::Register("GenBiasForceCollision");

  // -- build free flight operations:
  ConfigureForWorker();
}


void G4BOptrForceCollision::ConfigureForWorker()
{
  // -- start by remembering processes under biasing, and create needed biasing operations:
  if ( fSetup )
    {
      // -- get back ID created in master thread in case of MT (or reget just created ID in sequential mode):
      fForceCollisionModelID = G4PhysicsModelCatalog::Register("GenBiasForceCollision");

      const G4ProcessManager* processManager = fParticleToBias->GetProcessManager();
      const G4BiasingProcessSharedData* interfaceProcessSharedData = G4BiasingProcessInterface::GetSharedData( processManager );
      if ( interfaceProcessSharedData ) // -- sharedData tested, as is can happen a user attaches an operator
	                                // -- to a volume but without defining BiasingProcessInterface processes.
	{
	  for ( size_t i = 0 ; i < (interfaceProcessSharedData->GetPhysicsBiasingProcessInterfaces()).size(); i++ )
	    {
	      const G4BiasingProcessInterface* wrapperProcess =
		(interfaceProcessSharedData->GetPhysicsBiasingProcessInterfaces())[i];
	      G4String operationName = "FreeFlight-"+wrapperProcess->GetWrappedProcess()->GetProcessName();
	      fFreeFlightOperations[wrapperProcess] = new G4BOptnForceFreeFlight(operationName);
	    }
	}
      fSetup = false;
    }
}


void G4BOptrForceCollision::StartRun()
{
}


G4VBiasingOperation* G4BOptrForceCollision::ProposeOccurenceBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess)
{
  // -- does nothing if particle is not of requested type:
  if ( track->GetDefinition() != fParticleToBias ) return 0;

  // -- trying to get auxiliary track data...
  if ( fCurrentTrackData == nullptr )
    {
      // ... and if the track has no aux. track data, it means its biasing is not started yet (note that cloning is to happen first):
      fCurrentTrackData = (G4BOptrForceCollisionTrackData*)(track->GetAuxiliaryTrackInformation(fForceCollisionModelID));
      if ( fCurrentTrackData == nullptr ) return nullptr;
    }

  
  // -- Send force free flight to the callingProcess:
  // ------------------------------------------------
  // -- The track has been cloned in the previous step, it has now to be
  // -- forced for a free flight.
  // -- This track will fly with 0.0 weight during its forced flight:
  // -- it is to forbid the double counting with the force interaction track.
  // -- Its weight is restored at the end of its free flight, this weight
  // -- being its initial weight * the weight for the free flight travel,
  // -- this last one being per process. The initial weight is common, and is
  // -- arbitrary asked to the first operation to take care of it.
  if ( fCurrentTrackData->fForceCollisionState == ForceCollisionState::toBeFreeFlight )
    {
      G4BOptnForceFreeFlight* operation =  fFreeFlightOperations[callingProcess];
      if ( callingProcess->GetWrappedProcess()->GetCurrentInteractionLength() < DBL_MAX/10. )
	{
	  // -- the initial track weight will be restored only by the first DoIt free flight:
	  operation->ResetInitialTrackWeight(fInitialTrackWeight);
	  return operation;
	}
      else
	{
	  return nullptr;
	}
    }


  // -- Send force interaction operation to the callingProcess:
  // ----------------------------------------------------------
  // -- at this level, a copy of the track entering the volume was
  // -- generated (borned) earlier. This copy will make the forced
  // -- interaction in the volume.
  if ( fCurrentTrackData->fForceCollisionState == ForceCollisionState::toBeForced )
    {
      // -- Remember if this calling process is the first of the physics wrapper in
      // -- the PostStepGPIL loop (using default argument of method below):
      G4bool isFirstPhysGPIL = callingProcess-> GetIsFirstPostStepGPILInterface();
      
      // -- [*first process*] Initialize or update force interaction operation:
      if ( isFirstPhysGPIL )
	{
	  // -- first step of cloned track, initialize the forced interaction operation:
	  if ( track->GetCurrentStepNumber() == 1 ) fSharedForceInteractionOperation->Initialize( track );
	  else
	    {
	      if ( fSharedForceInteractionOperation->GetInitialMomentum() != track->GetMomentum() )
		{
		  // -- means that some other physics process, not under control of the forced interaction operation,
		  // -- has occured, need to re-initialize the operation as distance to boundary has changed.
		  // -- [ Note the re-initialization is only possible for a Markovian law. ]
		  fSharedForceInteractionOperation->Initialize( track );
		}
	      else
		{
		  // -- means that some other non-physics process (biasing or not, like step limit), has occured,
		  // -- but track conserves its momentum direction, only need to reduced the maximum distance for
		  // -- forced interaction.
		  // -- [ Note the update is only possible for a Markovian law. ]
		  fSharedForceInteractionOperation->UpdateForStep( track->GetStep() );
		}
	    }
	}
      
      // -- [*all processes*] Sanity check : it may happen in limit cases that distance to
      // -- out is zero, weight would be infinite in this case: abort forced interaction
      // -- and abandon biasing.
      if ( fSharedForceInteractionOperation->GetMaximumDistance() < DBL_MIN )
	{
	  fCurrentTrackData->Reset();
	  return 0;
	}
      
      // -- [* first process*] collect cross-sections and sample force law to determine interaction length
      // -- and winning process:
      if ( isFirstPhysGPIL )
	{
	  // -- collect cross-sections:
	  // -- ( Remember that the first of the G4BiasingProcessInterface triggers the update
	  // --   of these cross-sections )
	  const G4BiasingProcessSharedData* sharedData = callingProcess->GetSharedData();
	  for ( size_t i = 0 ; i < (sharedData->GetPhysicsBiasingProcessInterfaces()).size() ; i++ )
	    {
	      const G4BiasingProcessInterface* wrapperProcess = ( sharedData->GetPhysicsBiasingProcessInterfaces() )[i];
	      G4double interactionLength = wrapperProcess->GetWrappedProcess()->GetCurrentInteractionLength();
	      // -- keep only well defined cross-sections, other processes are ignored. These are not pathological cases
	      // -- but cases where a threhold effect par example (eg pair creation) may be at play. (**!**)
	      if ( interactionLength < DBL_MAX/10. )
		fSharedForceInteractionOperation->AddCrossSection( wrapperProcess->GetWrappedProcess(),  1.0/interactionLength );
	    }
	  // -- sample the shared law (interaction length, and winning process):
	  if ( fSharedForceInteractionOperation->GetNumberOfSharing() > 0 ) fSharedForceInteractionOperation->Sample();
	}
      
      // -- [*all processes*] Send operation for processes with well defined XS (see "**!**" above):
      G4VBiasingOperation* operationToReturn = nullptr;
      if ( callingProcess->GetWrappedProcess()->GetCurrentInteractionLength() < DBL_MAX/10. ) operationToReturn = fSharedForceInteractionOperation;
      return operationToReturn;


    } // -- end of "if ( fCurrentTrackData->fForceCollisionState == ForceCollisionState::toBeForced )"

  
  // -- other cases here: particle appearing in the volume by some
  // -- previous interaction : we decide to not bias these.
  return 0;
  
}


G4VBiasingOperation* G4BOptrForceCollision::ProposeNonPhysicsBiasingOperation(const G4Track* track,
									      const G4BiasingProcessInterface* /* callingProcess */)
{
  if ( track->GetDefinition() != fParticleToBias ) return nullptr;
  
  // -- Apply biasing scheme only to tracks entering the volume.
  // -- A "cloning" is done:
  // --  - the primary will be forced flight under a zero weight up to volume exit,
  // --    where the weight will be restored with proper weight for free flight
  // --  - the clone will be forced to interact in the volume.
  if ( track->GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ) // -- §§§ extent to case of a track shoot on the boundary ?
    {
      // -- check that track is free of undergoing biasing scheme ( no biasing data, or no active active )
      // -- Get possible track data:
      fCurrentTrackData = (G4BOptrForceCollisionTrackData*)(track->GetAuxiliaryTrackInformation(fForceCollisionModelID));
      if ( fCurrentTrackData != nullptr )
	{
	  if ( fCurrentTrackData->IsFreeFromBiasing() )
	    {
	      // -- takes "ownership" (some track data created before, left free, reuse it):
	      fCurrentTrackData->fForceCollisionOperator = this ;
	    }
	  else
	    {
	      // §§§ Would something be really wrong in this case ? Could this be that a process made a zero step ?
	    }
	}
      else
	{
	  fCurrentTrackData = new G4BOptrForceCollisionTrackData( this );
	  track->SetAuxiliaryTrackInformation(fForceCollisionModelID, fCurrentTrackData);
	}
      fCurrentTrackData->fForceCollisionState = ForceCollisionState::toBeCloned;
      fInitialTrackWeight = track->GetWeight();
      fCloningOperation->SetCloneWeights(0.0, fInitialTrackWeight);
      return fCloningOperation;
    }
  
  // -- 
  return nullptr;
}


G4VBiasingOperation* G4BOptrForceCollision::ProposeFinalStateBiasingOperation(const G4Track*, const G4BiasingProcessInterface* callingProcess)
{
  // -- Propose at final state generation the same operation which was proposed at GPIL level,
  // -- (which is either the force free flight or the force interaction operation).
  // -- count on the process interface to collect this operation:
  return callingProcess->GetCurrentOccurenceBiasingOperation();
}


void G4BOptrForceCollision::StartTracking( const G4Track* track )
{
  fCurrentTrack     = track;
  fCurrentTrackData = nullptr; 
}


void G4BOptrForceCollision::EndTracking()
{
  // -- check for consistency, operator should have cleaned the track:
  if ( fCurrentTrackData != nullptr )
    {
      if ( !fCurrentTrackData->IsFreeFromBiasing() )
	{
	  if ( (fCurrentTrack->GetTrackStatus() == fStopAndKill) || (fCurrentTrack->GetTrackStatus() == fKillTrackAndSecondaries) )
	    {
	      G4ExceptionDescription ed;
	      ed << "Current track deleted while under biasing by " << GetName() << ". Will result in inconsistencies.";
	      G4Exception(" G4BOptrForceCollision::EndTracking()",
			  "BIAS.GEN.18",
			  JustWarning,
			  ed);
	    }
	}
    } 
}


void G4BOptrForceCollision::OperationApplied( const G4BiasingProcessInterface*   callingProcess,
					      G4BiasingAppliedCase                          BAC,
					      G4VBiasingOperation*             operationApplied,
					      const G4VParticleChange*                          )
{
  
  if ( fCurrentTrackData == nullptr )
    {
      if ( BAC != BAC_None )
	{
	  G4ExceptionDescription ed;
	  ed << " Internal inconsistency : please submit bug report. " << G4endl;
	  G4Exception(" G4BOptrForceCollision::OperationApplied(...)",
		      "BIAS.GEN.20.1",
		      JustWarning,
		      ed); 
	}
      return;
    }
  
  if      ( fCurrentTrackData->fForceCollisionState == ForceCollisionState::toBeCloned )
    {
      fCurrentTrackData->fForceCollisionState = ForceCollisionState::toBeFreeFlight;
      auto cloneData                  = new G4BOptrForceCollisionTrackData( this );
      cloneData->fForceCollisionState = ForceCollisionState::toBeForced;
      fCloningOperation->GetCloneTrack()->SetAuxiliaryTrackInformation(fForceCollisionModelID, cloneData);
    }
  else if ( fCurrentTrackData->fForceCollisionState == ForceCollisionState::toBeFreeFlight )
    {
      if ( fFreeFlightOperations[callingProcess]->OperationComplete() ) fCurrentTrackData->Reset(); // -- off biasing for this track
    }
  else if ( fCurrentTrackData->fForceCollisionState == ForceCollisionState::toBeForced )
    {
      if ( operationApplied != fSharedForceInteractionOperation )
	{
	  G4ExceptionDescription ed;
	  ed << " Internal inconsistency : please submit bug report. " << G4endl;
	  G4Exception(" G4BOptrForceCollision::OperationApplied(...)",
		      "BIAS.GEN.20.2",
		      JustWarning,
		      ed); 
	}
      if ( fSharedForceInteractionOperation->GetInteractionOccured() )
	{
	  if ( operationApplied != fSharedForceInteractionOperation )
	    {
	      G4ExceptionDescription ed;
	      ed << " Internal inconsistency : please submit bug report. " << G4endl;
	      G4Exception(" G4BOptrForceCollision::OperationApplied(...)",
			  "BIAS.GEN.20.3",
			  JustWarning,
			  ed); 
	    } 
	}
    }
  else
    {
      if ( fCurrentTrackData->fForceCollisionState != ForceCollisionState::free )
	{
	  G4ExceptionDescription ed;
	  ed << " Internal inconsistency : please submit bug report. " << G4endl;
	  G4Exception(" G4BOptrForceCollision::OperationApplied(...)",
		      "BIAS.GEN.20.4",
		      JustWarning,
		      ed); 
	}
    }
}


void  G4BOptrForceCollision::OperationApplied( const G4BiasingProcessInterface*        /*callingProcess*/, G4BiasingAppliedCase                  /*biasingCase*/,
					       G4VBiasingOperation*         /*occurenceOperationApplied*/, G4double             /*weightForOccurenceInteraction*/,
					       G4VBiasingOperation*            finalStateOperationApplied, const G4VParticleChange*    /*particleChangeProduced*/ )
{
  
  if ( fCurrentTrackData->fForceCollisionState == ForceCollisionState::toBeForced )
    {
      if ( finalStateOperationApplied != fSharedForceInteractionOperation )
	{
	  G4ExceptionDescription ed;
	  ed << " Internal inconsistency : please submit bug report. " << G4endl;
	  G4Exception(" G4BOptrForceCollision::OperationApplied(...)",
		      "BIAS.GEN.20.5",
		      JustWarning,
		      ed); 
	}
      if ( fSharedForceInteractionOperation->GetInteractionOccured() ) fCurrentTrackData->Reset(); // -- off biasing for this track
    }
  else
    {
      G4ExceptionDescription ed;
      ed << " Internal inconsistency : please submit bug report. " << G4endl;
      G4Exception(" G4BOptrForceCollision::OperationApplied(...)",
		  "BIAS.GEN.20.6",
		  JustWarning,
		  ed);   
    }
  
}

