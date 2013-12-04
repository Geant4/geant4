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
#include "G4BiasingProcessInterface.hh"

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

G4BOptrForceCollision::G4BOptrForceCollision(G4String particleName, G4String name)
  : G4VBiasingOperator(name),
    fFirstProcess(0),   fLastProcess(0),
    fSetup(true)
{
  fSharedForceInteractionOperation = new G4BOptnForceCommonTruncatedExp("SharedForceInteraction");
  fCloningOperation                = new G4BOptnCloning("Cloning");
  fParticle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  
  if ( fParticle == 0 )
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
    fFirstProcess(0),   fLastProcess(0),
    fSetup(true)
{
  fSharedForceInteractionOperation = new G4BOptnForceCommonTruncatedExp("SharedForceInteraction");
  fCloningOperation                = new G4BOptnCloning("Cloning");
  fParticle                        = particle;
}

G4BOptrForceCollision::~G4BOptrForceCollision()
{
  for ( std::map< const G4BiasingProcessInterface*, G4BOptnForceFreeFlight* >::iterator it = fFreeFlightOperations.begin() ;
	it != fFreeFlightOperations.end() ;
	it++ ) delete (*it).second;
  delete fSharedForceInteractionOperation;
  delete fCloningOperation;
}


G4VBiasingOperation* G4BOptrForceCollision::ProposeOccurenceBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess)
{
  if ( track->GetDefinition() != fParticle ) return 0;

  // -- start by remembering processes under biasing, and create needed biasing operations:
  // -- ( Might consider moving this in a less called method. )
  if ( fSetup )
    {
      if ( ( fFirstProcess == 0 ) && ( callingProcess->GetIsFirstPostStepGPILInterface() ) ) fFirstProcess = callingProcess;
      if ( fLastProcess == 0 )
	{
	  fProcesses.push_back(callingProcess);
	  G4String operationName = "FreeFlight-"+callingProcess->GetWrappedProcess()->GetProcessName();
	  fFreeFlightOperations[callingProcess] = new G4BOptnForceFreeFlight(operationName);
	  if ( callingProcess->GetIsLastPostStepGPILInterface() )
	    {
	      fLastProcess = callingProcess;
	      fSetup = false;
	    }
	}
    }
  
  
  // -- Send force interaction operation to the callingProcess:
  // ----------------------------------------------------------
  // -- at this level, a copy of the track entering the volume was
  // -- generated (borned) earlier. This copy will make the forced
  // -- interaction in the volume.
  if ( GetBirthOperation( track ) == fCloningOperation )
    {
      // -- forced interaction already occured, abort.
      // -- [Note that if ones would redo a forced interaction, it
      // -- would require an other cloning before, if one wants to conserve
      // -- the weight.]
      if ( fSharedForceInteractionOperation->GetInteractionOccured() ) return 0;
      
      if ( callingProcess == fFirstProcess )
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
      
      // -- Sanity check : it may happen in limit cases that distance to out is zero,
      // -- weight would be infinite in this case. Abort forced interaction.
      if ( fSharedForceInteractionOperation->GetMaximumDistance() < DBL_MIN ) return 0;
      
      // -- conditions to apply forced interaction are met, set up physics:
      // -- Collects well-defined cross-sections...
      G4double      currentInteractionLength =  callingProcess->GetWrappedProcess()->GetCurrentInteractionLength();
      G4VBiasingOperation* operationToReturn = 0;
      if ( currentInteractionLength < DBL_MAX/10. )
	{
	  fSharedForceInteractionOperation->AddCrossSection( callingProcess->GetWrappedProcess(), 1.0/currentInteractionLength );
	  operationToReturn = fSharedForceInteractionOperation;
	}
      // -- ... if current process is the last one, cross-section collection is finished.
      // -- Operation is ready to sample the force interaction law and to randomly select
      // -- the process to be applied.
      if ( ( callingProcess == fLastProcess ) &&
	   ( fSharedForceInteractionOperation->GetCommonTruncatedExpLaw()->GetNumberOfSharing() > 0 ) )
	fSharedForceInteractionOperation->Sample();
      
      // -- back to current process (last one or not), if meaningful cross section, bias it returning the force operation:
      return operationToReturn;

    } // -- end of " if ( GetBirthOperation( track ) == fCloningOperation ) "
  
  
  
  // -- Send force free flight to the callingProcess:
  // ------------------------------------------------
  // -- At this point a track cloning happened in the previous step,
  // -- we request the continuing/current track to be forced flight.
  // -- Note this track will fly with 0.0 weight during its forced flight:
  // -- it is to forbid the double counting with the force interaction track.
  // -- Its weight is restored at the end of its free flight, this weight
  // -- being its initial weight * the weight for the free flight travel,
  // -- this last one being per process. The initial weight is common, and is
  // -- arbitrary asked to the first operation to take care of it.
  if ( fPreviousOperationApplied == fCloningOperation )
    {
      G4BOptnForceFreeFlight* operation =  fFreeFlightOperations[callingProcess];
      if ( callingProcess->GetWrappedProcess()->GetCurrentInteractionLength() < DBL_MAX/10. )
	{
	  // -- the initial track weight will be restored only by the first DoIt free flight:
	  operation->ResetInitialTrackWeight(fInitialTrackWeight);
	  return operation;
	}
    }
  // -- If forced flight was already requested for the calling process, this flight
  // -- may have been interupted by some other non-physics process : request
  // -- process to continue with same forced flight:
  // -- ** ?? ** Incorrect here in case of flight interrupted by a *physics process*. Ie case
  // -- ** ?? ** of a subset of physics processes forced is not treated correctly here ** ?? **
  else if ( callingProcess->GetPreviousOccurenceBiasingOperation() == fFreeFlightOperations[callingProcess] )
    return fFreeFlightOperations[callingProcess];
  
  // -- other cases here: particle appearing in the volume by some
  // -- previous interaction : we decide to not bias these.
  return 0;
  
}


G4VBiasingOperation* G4BOptrForceCollision::ProposeNonPhysicsBiasingOperation(const G4Track* track,
									      const G4BiasingProcessInterface*)
{
  if ( track->GetDefinition() != fParticle ) return 0;
  
  // -- bias only tracks entering the volume.
  // -- A "cloning" is done:
  // --  - the primary will be forced flight under a zero weight up to volume exit,
  // --    where the weight will be restored with proper weight for free flight
  // --  - the clone will be forced to interact in the volume.
  if ( track->GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary )
    {
      fSharedForceInteractionOperation->SetInteractionOccured( false );
      fInitialTrackWeight = track->GetWeight();
      fCloningOperation->SetCloneWeights(0.0, fInitialTrackWeight);
      return fCloningOperation;
    }

  // -- for all other cases: does nothing
  return 0;
}


void G4BOptrForceCollision::StartTracking( const G4Track* )
{
  fPreviousOperationApplied = 0;
}


void G4BOptrForceCollision::ExitBiasing( const G4Track* track, const G4BiasingProcessInterface* )
{
  fPreviousOperationApplied = 0;
  // -- clean-up track for what happens with this operator:
  ForgetTrack ( track );
}


void G4BOptrForceCollision::OperationApplied( const G4BiasingProcessInterface*   callingProcess, G4BiasingAppliedCase,
					      G4VBiasingOperation*             operationApplied, const G4VParticleChange* particleChangeProduced )
{
  fPreviousOperationApplied = operationApplied;
  if ( operationApplied == fCloningOperation )
    RememberSecondaries( callingProcess, operationApplied, particleChangeProduced );
}
