// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastSimulationManagerProcess.cc,v 1.3 1999-04-28 10:06:43 mora Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4FastSimulationProcess.cc
//
//  Description:
//    The process that triggers the parameterised simulations,
//    if any.
//
//  History:
//    August 97: First implementation. Verderi && MoraDeFreitas.
//---------------------------------------------------------------

#include "G4ios.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleChange.hh"

//----------------------------
// Constructor with only name:
//----------------------------
G4FastSimulationManagerProcess::
G4FastSimulationManagerProcess(const G4String& processName,
			       G4ProcessType theType) : 
  G4VProcess(processName,theType)
{
  pParticleChange = &aDummyParticleChange;

  fGhostTouchable       = new G4TouchableHistory();
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created " << endl;
  }
  fGhostFieldPropagator = new G4PropagatorInField(&fGhostNavigator,
						   G4TransportationManager::
						   GetTransportationManager()->
						   GetFieldManager());
}

// -----------
// Destructor:
// -----------
G4FastSimulationManagerProcess::~G4FastSimulationManagerProcess()
{
  delete fGhostTouchable;
  fGhostTouchable = 0;
  delete fGhostFieldPropagator;
  fGhostFieldPropagator = 0;
}


//----------------------------------------------------------
//
//  PostStepGetPhysicalInteractionLength()
//
//    This method is used to trigger the parameterised 
//    simulation.
//
//----------------------------------------------------------
//
//    To be triggered the conditions that must be
//    fullfilled are :
//
//    1) a call to track->GetVolume()->GetLogicalVolume()->
//    GetFastSimulationManager() returns a G4FastSimulationManager
//    not NULL pointer (kept in fFastSimulationManager);
//
//    AND
//
//    2)a call to this fFastSimulationManager->GetTrigger() method 
///     returns a G4bool True. 
//    
//
//    If the fFastSimulationManager* is triggered then:
//        
//            * condition = ExclusivelyForced 
//            * returns 0.0
//
//    else:
//            * condition = NotForced
//            * returns DBL_MAX
//
//-----------------------------------------------------------
G4double 
G4FastSimulationManagerProcess::PostStepGetPhysicalInteractionLength(
         const G4Track& track, 
         G4double   previousStepSize, 
         G4ForceCondition* condition)
{
  // ------------------------------
  // Prelude for initial Step only:
  // ------------------------------
  // ----------------------------------------------------------
  // -- If starting a track, check if there is a parallel World
  // -- for the new particle being tracked.
  // -- (may be Done more than once in very special cases,
  // -- since the fStarTracking is switched off lower.)
  // ----------------------------------------------------------
  if(fStartTracking)
    {
      G4VFlavoredParallelWorld* flavoredWorld = 
	G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->
	GetFlavoredWorldForThis(track.GetDynamicParticle()->GetDefinition());

      fGhostWorld               = 0;
      if (flavoredWorld)
	fGhostWorld             = flavoredWorld->GetThePhysicalVolumeWorld();

      if (fGhostWorld)
	{
	  fGhostNavigator.SetWorldVolume(fGhostWorld);
	  fUpdateGhostTouchable = true;
	  fGhostSafety          = -1.;
	  // Navigation in Field:
	  fParticleCharge       = track.GetDynamicParticle()->GetDefinition()->GetPDGCharge();
	}
      else
	{
	  fOutOfGhostWorld      = true;
	}
    }

  // -------------------------------------------------------
  //
  //     Parameterisation Trigger invokation sequence :
  //
  // -------------------------------------------------------
  fFastSimulationTrigger = false;

  //---------------------------------------------
  // Normal Dispatcher (for tracking geometry) :
  //---------------------------------------------
  if(fFastSimulationManager = 
     track.GetVolume()->GetLogicalVolume()->GetFastSimulationManager())
    {
      // Yes, so should us trigger a fast simulation model now ?
      if(fFastSimulationTrigger = 
	 fFastSimulationManager->PostStepGetFastSimulationManagerTrigger(track))
	{
	  // Yes, lets take the control over the stepping !
	  *condition = ExclusivelyForced;
	  return 0.0;
	}
    }
  
  //-----------------------------------------
  // * Parallel geometry Dispatcher if any.
  // * If no trigger in the parallel geometry
  // returns the Conditionally forced signal.
  //-----------------------------------------
  if(fGhostWorld)
    {
      // First, update the touchable history of ghost if necessary, ie when:
      //      - this is first Step of tracking
      //      - the last Step was limited by ghost geometry
      // otherwise performs a single Locate
      if (fUpdateGhostTouchable)
	{
	  fUpdateGhostTouchable=false;
	  if (fStartTracking)
	    {
	      fGhostNavigator.LocateGlobalPointAndUpdateTouchable(
					          track.GetPosition(),
						  track.GetMomentumDirection(),
					          fGhostTouchable,
						  false);
	      fStartTracking = false;
	    }
	  else
	    {
 	      fGhostNavigator.SetGeometricallyLimitedStep();
	      fGhostNavigator.LocateGlobalPointAndUpdateTouchable(
						  track.GetPosition(),
						  track.GetMomentumDirection(),
						  fGhostTouchable);
	    }
	  fOutOfGhostWorld = (fGhostTouchable->GetVolume() == 0);
	}
      else if (previousStepSize > 0.0)
	{
	  // G4ThreeVector direction= track.GetMomentumDirection();
	  // fGhostNavigator.LocateGlobalPointAndSetup(track.GetPosition(), &direction, true);
	  // The above looks not enough in case of mag-field, not clear why ?!?
	  fGhostNavigator.LocateGlobalPointAndUpdateTouchable(
					      track.GetPosition(),
					      track.GetMomentumDirection(),
					      fGhostTouchable);
	}


      if (!fOutOfGhostWorld)
	{
	  if(fFastSimulationManager = 
	     fGhostTouchable->GetVolume()->GetLogicalVolume()->GetFastSimulationManager())
	    {
	      // Yes, so should us trigger a fast simulation model now ?
	      if(fFastSimulationTrigger = 
		 fFastSimulationManager->PostStepGetFastSimulationManagerTrigger(track,&fGhostNavigator))
		{
		  // Yes, lets take the control over the stepping !
		  *condition = ExclusivelyForced;
		  return 0.0;
		} 
	    }
	  // No ghost trigger has occured (ie no manager or a manager but no trigger)
	  // Informs the stepping that PostStepDoIt may be called if the AlongGPIL
	  // will limit the Step.
	  // PostStepDoIt will switch on fUpdateGhostTouchable flag,
	  // and eventually correct position and mometum in case of
	  // mag-field.
	  *condition = Conditionally;
	  return DBL_MAX;
	}
    }
  
  *condition = NotForced;
  return DBL_MAX;
}

//------------------------------------
//
//             PostStepDoIt()
//
//------------------------------------
G4VParticleChange* G4FastSimulationManagerProcess::PostStepDoIt(
     const G4Track& track,
     const G4Step&  Step)
{
  if (fFastSimulationTrigger) 
    {
      // Executes the ParameterisedSimulation code.
      // and gets the ParameterisedSimulation response.
      G4VParticleChange* Response=
	fFastSimulationManager->InvokePostStepDoIt();

      // If the particle is still alive, suspend it
      // to re-initialise the other process.
      if (Response->GetStatusChange() != fStopAndKill)
	Response->SetStatusChange(fSuspend);
      // Returns the Response
      return Response;
    }
  else
    {
      if (fOutOfGhostWorld)   G4cout << " ????????????? Problem of logic in GHOST param !!! " << endl;
      // No FastSimulationManager has asked for trigger. That means here:
      //       - the PostStep is "Conditionally" forced because
      //       - the Step has been limited by the ALong of G4FSMP because
      //       - the track has reached a ghost boundary
      // => the touchable history must be updated.
      // fUpdateGhostTouchable = true; : MOVED BELOW: COMPLICATIONS WITH MAG-FIELD !!

      if (fFieldExertsForce)
	{
	  // -----------------------------------------------------------
	  // in case of field, correct for small discrepancy on position
	  // and momentum computed by the Transportation.
	  // HOWEVER, this correction can be somewhat problematic when 
	  // the end point, which is close here to a Ghost boundary by
	  // construction, is also close to a boundary of the tracking
	  // geometry.
	  // Thus we correct when the safety in the tracking geometry
	  // allows to move the point enough. Otherwise we issue
	  // a warning.
	  // ----------------------------------------------------------
	  fTrackingNavigator.SetWorldVolume(G4TransportationManager::GetTransportationManager()->
					    GetNavigatorForTracking()->GetWorldVolume());
	  fTrackingNavigator.LocateGlobalPointAndUpdateTouchable(
						  track.GetPosition(),	  
						  track.GetMomentumDirection(),
					          &fTrackingHistory);
	  G4VPhysicalVolume* trackingVolume = fTrackingHistory.GetVolume();

	  G4double trackingSafety(0.0);
	  G4double trackingLinearDistance = 
	    fTrackingNavigator.ComputeStep(
			   track.GetPosition(),
			   track.GetMomentumDirection(),
			   DBL_MAX,
			   trackingSafety);
	  
	  xParticleChange.Initialize(track);
	  G4ThreeVector deltaPos = fGhostFieldPropagator->EndPosition() - track.GetPosition();
	  if (trackingSafety < deltaPos.mag())
	    {
	      fUpdateGhostTouchable = true; // What should we do here ?!?
	                                    // false => lot a microscopic steps
	                                    // is true rasonnably safe ?
	      G4cout << GetProcessName() << " : WARNING: Can not correct for\ndifference between tracking and ghost mag-field computations." << endl;
	    }
	  else
	    {
	      // easy case where we have enough room to displace the point where we need:
	      fUpdateGhostTouchable = true;
	      xParticleChange.SetPositionChange(fGhostFieldPropagator->EndPosition());
	      xParticleChange.SetMomentumChange(fGhostFieldPropagator->EndMomentumDir());
	    }
	  return &xParticleChange;
	}
      else
	{
	  fUpdateGhostTouchable = true;
	  pParticleChange->Initialize(track);
	  return pParticleChange; 
	}
    }
}


G4double G4FastSimulationManagerProcess::AlongStepGetPhysicalInteractionLength(
				              const G4Track& track,
				              G4double  previousStepSize,
				              G4double  currentMinimumStep,
				              G4double& proposedSafety,
				              G4GPILSelection* selection
				              )
{
  // Informs the stepping that G4FastsimulationProcess wants
  // to be able to limit the Step.
  *selection = CandidateForSelection;

  G4double returnedStep = DBL_MAX;

  if (!fOutOfGhostWorld)
    {
      if (previousStepSize > 0.) fGhostSafety -= previousStepSize;
      else                       fGhostSafety = -1.;
      if (fGhostSafety < 0.) fGhostSafety = 0.0;
      
      // ------------------------------------------
      // Determination of the proposed STEP LENGTH:
      // ------------------------------------------
      if (currentMinimumStep <= fGhostSafety)
	{
	  returnedStep = currentMinimumStep;
	}
      else // (currentMinimumStep > fGhostSafety: GFSMP may limit the Step)
	{
	  fFieldExertsForce = (fParticleCharge != 0.0) &&
	                      (G4TransportationManager::GetTransportationManager()->
	                       GetFieldManager()->DoesFieldExist());
	  if ( !fFieldExertsForce )
	    {
	      fGhostStepLength = fGhostNavigator.ComputeStep(
						 track.GetPosition(),
						 track.GetMomentumDirection(),
						 currentMinimumStep,
						 fGhostSafety);

	    }
	  else
	    {
	      fGhostFieldPropagator->SetChargeMomentumMass(
						 fParticleCharge,
						 track.GetDynamicParticle()->
						 GetTotalMomentum(),
						 track.GetDynamicParticle()->
						 GetMass());

	      fGhostStepLength = fGhostFieldPropagator->ComputeStep( 
						        track.GetPosition(), 
						        track.GetMomentumDirection(),
						        currentMinimumStep,
						        fGhostSafety);
	      }
	  fPreSafety           = fGhostSafety;
	  returnedStep         = fGhostStepLength;
	}

      // ----------------------------------------------
      // Returns the fGhostSafety as the proposedSafety
      // The SteppingManager will take care of keeping
      // the smallest one.
      // ----------------------------------------------
      proposedSafety           = fGhostSafety;
    }
  return returnedStep;
}

G4VParticleChange* G4FastSimulationManagerProcess::AlongStepDoIt(
				             const G4Track& track,
				             const G4Step&  Step
				             )
{
  // Dummy ParticleChange ie: does nothing
  pParticleChange->Initialize(track);
  return pParticleChange; 
}

//--------------------------------------------
//
//         At Rest parameterisation:
//
//--------------------------------------------
//
//   AtRestGetPhysiscalInteractionLength:
//
//--------------------------------------------
G4double 
G4FastSimulationManagerProcess::AtRestGetPhysicalInteractionLength(
         const G4Track& track, 
         G4ForceCondition* condition)
{
  if ( (fFastSimulationManager =
	track.GetVolume()->GetLogicalVolume()->GetFastSimulationManager())
       !=  0 )
    if (fFastSimulationManager->AtRestGetFastSimulationManagerTrigger(track))
      {
	// "*condition = ExclusivelyForced;" Not yet available
	// for AtRest actions, so we use the trick below. However
	// it is not garantee to work in the situation the track is
	// alive after parameterisation.
	// NOTE: Problem: the present stepping doesn't care if the particle
	// has been killed between two AtRestDoIt() !!!
	*condition = NotForced;
	return -1.0; // TEMPORARY TRICK TO TAKE CONTROL !!!!
      }
  
  //-----------------------------------------
  // * Parallel geometry Dispatcher if any.
  // * If no trigger in the parallel geometry
  // * returns DBL_MAX and NotForced signal.
  //-----------------------------------------
  if(fGhostWorld)
    {
      // First, update the touchable history of ghost if necessary, ie when:
      //      - this is first Step of tracking
      //      - the last Step was limited by ghost geometry
      // otherwise performs a single Locate
      if (fUpdateGhostTouchable)
	{
	  fUpdateGhostTouchable=false;
	  if (fStartTracking)
	    {
	      fGhostNavigator.LocateGlobalPointAndUpdateTouchable(
						  track.GetPosition(),
					          fGhostTouchable,
						  false);
	      fStartTracking = false;
	    }
	  else
	    {
	      fGhostNavigator.SetGeometricallyLimitedStep();
	      fGhostNavigator.LocateGlobalPointAndUpdateTouchable(
						  track.GetPosition(),	  
					          fGhostTouchable);

	    }
	  fOutOfGhostWorld = (fGhostTouchable->GetVolume() == 0);
	}

      if (!fOutOfGhostWorld)
	{
	  if(fFastSimulationManager = 
	     fGhostTouchable->GetVolume()->GetLogicalVolume()->GetFastSimulationManager())
	    {
	      // Should it trigger a fast simulation model now ?
	      if(fFastSimulationTrigger = 
		 fFastSimulationManager->
		 AtRestGetFastSimulationManagerTrigger(track,&fGhostNavigator))
		{
		  // "*condition = ExclusivelyForced;" Not yet available
		  // for AtRest actions, so we use the trick below. However
		  // it is not garantee to work in the situation the track is
		  // alive after parameterisation.
		  // NOTE: Problem: the present stepping doesn't care if the particle
		  // has been killed between two AtRestDoIt() !!!
		  *condition = NotForced;
		  return -1.0; // TEMPORARY TRICK TO TAKE CONTROL !!!!
		} 
	    }
	}
    }

  // No, avoid to take control over the PostStepDoIt methods.
  *condition = NotForced;
  return DBL_MAX;
}

//-----------------------------------------------
//
//                  AtRestDoIt:
//
//-----------------------------------------------
G4VParticleChange* G4FastSimulationManagerProcess::AtRestDoIt(
     const G4Track& track,
     const G4Step& Step)
{
  return fFastSimulationManager->InvokeAtRestDoIt();
}

void G4FastSimulationManagerProcess::StartTracking() 
{
  fStartTracking = true;
}

void G4FastSimulationManagerProcess::Verbose() const
{
  /*  G4cout << "     >>>>> Trigger Status : ";
  switch(fFastSimulationManager->GetTriggerStatus())
    {
    case NoModel:
      G4cout << "NoModel" << endl;
      break;
    case OnBoundaryButLeaving:
      G4cout << "OnBoundaryButLeaving" << endl;
      break;
    case OneModelTrigger:
      G4cout << "OneModelTrigger" << endl;
      break;
    case NoModelTrigger:
      G4cout << "NoModelTrigger" << endl;
      break;
    case Undefined:
      G4cout << "Undefined" << endl;
      break;
    default:
      G4cout << " Bizarre..." << endl;
      break;
    }*/
}
