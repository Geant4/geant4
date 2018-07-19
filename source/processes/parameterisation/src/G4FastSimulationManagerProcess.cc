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
//
// $Id: G4FastSimulationManagerProcess.cc 101152 2016-11-08 08:07:39Z gcosmo $
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
//    October 06: move to parallel geometry scheme, M. Verderi
//---------------------------------------------------------------

#include "G4ios.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4TransportationManager.hh"
#include "G4PathFinder.hh"
#include "G4ParticleChange.hh"
#include "G4FieldTrackUpdator.hh"

#define  PARANOIA

G4FastSimulationManagerProcess::
G4FastSimulationManagerProcess(const G4String& processName,
			       G4ProcessType       theType) : 
  G4VProcess(processName,theType),
  fWorldVolume          ( nullptr ),
  fIsTrackingTime       ( false   ),
  fIsFirstStep          ( false   ),
  fGhostNavigator       ( nullptr ),
  fGhostNavigatorIndex  ( -1      ),
  fIsGhostGeometry      ( false   ),
  fGhostSafety          ( -1.0    ),
  fFieldTrack           ( '0'     ),
  fFastSimulationManager( nullptr ),
  fFastSimulationTrigger( false   )
{
  // -- set Process Sub Type
  SetProcessSubType(static_cast<int>(FASTSIM_ManagerProcess));


  fPathFinder            = G4PathFinder::GetInstance();
  fTransportationManager = G4TransportationManager::GetTransportationManager();
  
  SetWorldVolume(fTransportationManager->GetNavigatorForTracking()->GetWorldVolume()->GetName());
  if (verboseLevel>0) G4cout << "G4FastSimulationManagerProcess `" << GetProcessName() 
			     << "' is created, and will message geometry with world volume `" 
			     << fWorldVolume->GetName() << "'." << G4endl;
  G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->AddFSMP(this);
}


G4FastSimulationManagerProcess::
G4FastSimulationManagerProcess(const G4String&     processName,
			       const G4String& worldVolumeName,
			       G4ProcessType           theType) :
  G4VProcess(processName,theType),
  fWorldVolume          ( nullptr ),
  fIsTrackingTime       ( false   ),
  fIsFirstStep          ( false   ),
  fGhostNavigator       ( nullptr ),
  fGhostNavigatorIndex  ( -1      ),
  fIsGhostGeometry      ( false   ),
  fGhostSafety          ( -1.0    ),
  fFieldTrack           ( '0'     ),
  fFastSimulationManager( nullptr ),
  fFastSimulationTrigger( false   )
{
  // -- set Process Sub Type
  SetProcessSubType(static_cast<int>(FASTSIM_ManagerProcess));


  fPathFinder            = G4PathFinder::GetInstance();
  fTransportationManager = G4TransportationManager::GetTransportationManager();

  SetWorldVolume(worldVolumeName);
  if (verboseLevel>0) G4cout << "G4FastSimulationManagerProcess `" << GetProcessName() 
			     << "' is created, and will message geometry with world volume `" 
			     << fWorldVolume->GetName() << "'." << G4endl;
  G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->AddFSMP(this);
}


G4FastSimulationManagerProcess::
G4FastSimulationManagerProcess(const G4String&    processName,
			       G4VPhysicalVolume* worldVolume,
			       G4ProcessType      theType) :
  G4VProcess(processName,theType),
  fWorldVolume          ( nullptr ),
  fIsTrackingTime       ( false   ),
  fIsFirstStep          ( false   ),
  fGhostNavigator       ( nullptr ),
  fGhostNavigatorIndex  ( -1      ),
  fIsGhostGeometry      ( false   ),
  fGhostSafety          ( -1.0    ),
  fFieldTrack           ( '0'     ),
  fFastSimulationManager( nullptr ),
  fFastSimulationTrigger( false   )
{
  // -- set Process Sub Type
  SetProcessSubType(static_cast<int>(FASTSIM_ManagerProcess));
  

  fPathFinder            = G4PathFinder::GetInstance();
  fTransportationManager = G4TransportationManager::GetTransportationManager();

  SetWorldVolume(worldVolume);
  if (verboseLevel>0) G4cout << "G4FastSimulationManagerProcess `" << GetProcessName() 
			     << "' is created, and will message geometry with world volume `" 
			     << fWorldVolume->GetName() << "'." << G4endl;
  G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->AddFSMP(this);
}


G4FastSimulationManagerProcess::~G4FastSimulationManagerProcess()
{
  G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->RemoveFSMP(this);
}


// -----------------------
//   User access methods:
// -----------------------
void G4FastSimulationManagerProcess::SetWorldVolume(G4String newWorldName)
{
  if (fIsTrackingTime)
    {
      G4ExceptionDescription ed;
      ed << "G4FastSimulationManagerProcess `" << GetProcessName()
	 << "': changing of world volume at tracking time is not allowed." << G4endl;
      G4Exception("G4FastSimulationManagerProcess::SetWorldVolume(const G4String)",
		  "FastSim002",
		  JustWarning, ed,
		  "Call ignored.");
    }
  else
    {
      G4VPhysicalVolume* newWorld = fTransportationManager->IsWorldExisting(newWorldName);
      if (newWorld == 0)
	{
	  G4ExceptionDescription  tellWhatIsWrong;
	  tellWhatIsWrong << "Volume newWorldName = `" <<  newWorldName
			  << "' is not a parallel world nor the mass world volume."
			  << G4endl;
	  G4Exception("G4FastSimulationManagerProcess::SetWorldVolume(const G4String)",
		      "FastSim003",
		      FatalException,
		      tellWhatIsWrong);
	}
      if (verboseLevel>0)
      {
	if (fWorldVolume) G4cout << "G4FastSimulationManagerProcess `" << GetProcessName()
				 << "': changing world volume from '"  << fWorldVolume->GetName() 
				 << "' to `" << newWorld << "'." << G4endl;
	else              G4cout << "G4FastSimulationManagerProcess `" << GetProcessName()
				 << "': setting world volume from to `"<< newWorld->GetName() << "'." << G4endl;
      }
      fWorldVolume = newWorld;
    }
}


void G4FastSimulationManagerProcess::SetWorldVolume(G4VPhysicalVolume* newWorld)
{
  if (newWorld) SetWorldVolume(newWorld->GetName());
  else
    {
      G4ExceptionDescription  tellWhatIsWrong;
      tellWhatIsWrong << "Null pointer passed for world volume." << G4endl;
      G4Exception("G4FastSimulationManagerProcess::SetWorldVolume(const G4VPhysicalVolume* newWorld)",
		  "FastSim004",
		  FatalException,
		  tellWhatIsWrong); 
    }
}


// --------------------
//  Start/End tracking:
// --------------------
void G4FastSimulationManagerProcess::StartTracking(G4Track* track)
{
  fIsTrackingTime = true;
  fIsFirstStep    = true;
  
  // -- fetch the navigator (and its index) and activate it:
  G4TransportationManager* transportationManager = G4TransportationManager::GetTransportationManager();
  fGhostNavigator                            = transportationManager->GetNavigator(fWorldVolume);
  fIsGhostGeometry                           = (fGhostNavigator != transportationManager->GetNavigatorForTracking());
  if (fIsGhostGeometry) fGhostNavigatorIndex = transportationManager->ActivateNavigator(fGhostNavigator);
  else                  fGhostNavigatorIndex = -1;

  fPathFinder->PrepareNewTrack(track->GetPosition(), track->GetMomentumDirection());
}


void
G4FastSimulationManagerProcess::
EndTracking()
{
  fIsTrackingTime = false;
  if ( fIsGhostGeometry ) fTransportationManager->DeActivateNavigator(fGhostNavigator); 
}


// ------------------------------------------
//   PostStepGetPhysicalInteractionLength():
// ------------------------------------------
G4double 
G4FastSimulationManagerProcess::
PostStepGetPhysicalInteractionLength(const G4Track&               track, 
				     G4double,
				     G4ForceCondition*        condition)
{
  // -- Get current volume, and check for presence of fast simulation manager.
  // -- For the case of the navigator for tracking (fGhostNavigatorIndex == 0)
  // -- we use the track volume. This allows the code to be valid for both
  // -- cases where the PathFinder is used (G4CoupledTranportation) or not
  // -- (G4Transportation).
  const G4VPhysicalVolume* currentVolume(0);
  if ( fIsGhostGeometry )  currentVolume = fPathFinder->GetLocatedVolume(fGhostNavigatorIndex);
  else                     currentVolume = track.GetVolume();

  if ( currentVolume )
    {
      fFastSimulationManager = currentVolume->GetLogicalVolume()->GetFastSimulationManager();
      if( fFastSimulationManager )
	{
	  // Ask for trigger:
	  fFastSimulationTrigger = fFastSimulationManager->PostStepGetFastSimulationManagerTrigger(track, fGhostNavigator);
	  if( fFastSimulationTrigger )
	    {
	      // Take control over stepping:
	      *condition = ExclusivelyForced;
	      return 0.0;
	    }
	}     
    }
  
  // -- no fast simulation occuring there:
  *condition = NotForced;
  return DBL_MAX;
}

//------------------------------------
//             PostStepDoIt()
//------------------------------------
G4VParticleChange*
G4FastSimulationManagerProcess::
PostStepDoIt(const G4Track&,
	     const G4Step&)
{
  G4VParticleChange* finalState = fFastSimulationManager->InvokePostStepDoIt();
  
  // If the particle is still alive, suspend it to force physics re-initialisation:
  if (finalState->GetTrackStatus() != fStopAndKill) finalState->ProposeTrackStatus(fSuspend);
  
  return finalState;
}


G4double 
G4FastSimulationManagerProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&                track,
				      G4double           previousStepSize,
				      G4double         currentMinimumStep, 
				      G4double&            proposedSafety, 
				      G4GPILSelection*          selection)
{
  
  *selection            = NotCandidateForSelection;
  G4double returnedStep = DBL_MAX;

  // ---------------------------------------------------
  // -- Below code valid for ghost geometry, otherwise
  // -- useless for fast simulation attached to mass
  // -- geometry. Warn user in case along used for 
  // -- mass geometry ?
  // --------------------------------------------------
  if ( fIsGhostGeometry )
    {
      static G4ThreadLocal G4FieldTrack *endTrack_G4MT_TLS_ = 0 ;
      if (!endTrack_G4MT_TLS_) endTrack_G4MT_TLS_ = new  G4FieldTrack ('0') ;
      G4FieldTrack &endTrack = *endTrack_G4MT_TLS_;
      
      static G4ThreadLocal ELimited *eLimited_G4MT_TLS_ = 0 ;
      if (!eLimited_G4MT_TLS_) eLimited_G4MT_TLS_ = new  ELimited  ;
      ELimited &eLimited = *eLimited_G4MT_TLS_;
      
      if (previousStepSize > 0.) fGhostSafety -= previousStepSize;
      if (fGhostSafety < 0.)     fGhostSafety = 0.0;
      
      // ------------------------------------------
      // Determination of the proposed step length:
      // ------------------------------------------
      if (currentMinimumStep <= fGhostSafety && currentMinimumStep > 0.)
	{
	  // -- No chance to limit the step, as proposed move inside safety
	  returnedStep   = currentMinimumStep;
	  proposedSafety = fGhostSafety - currentMinimumStep;
	}
      else
	{
	  // -- Proposed move exceeds safety, need to state
	  G4FieldTrackUpdator::Update(&fFieldTrack, &track);
	  returnedStep = fPathFinder->ComputeStep(fFieldTrack,
						  currentMinimumStep,
						  fGhostNavigatorIndex,
						  track.GetCurrentStepNumber(),
						  fGhostSafety,
						  eLimited,
						  endTrack,
						  track.GetVolume());
	  
	  if(eLimited == kDoNot) fGhostSafety = fGhostNavigator->ComputeSafety(endTrack.GetPosition()); // -- step no limited by ghost
	  proposedSafety = fGhostSafety;
	  if      (eLimited == kUnique || eLimited == kSharedOther) *selection     = CandidateForSelection;
	  else if (eLimited == kSharedTransport)                     returnedStep *= (1.0 + 1.0e-9);   // -- Expand to disable its selection in Step Manager comparison
	}
    }
  

  // ----------------------------------------------
  // Returns the fGhostSafety as the proposedSafety
  // The SteppingManager will take care of keeping
  // the smallest one.
  // ----------------------------------------------
  return returnedStep;
}

G4VParticleChange*
G4FastSimulationManagerProcess::
AlongStepDoIt(const G4Track& track,
	      const G4Step&)
{
  fDummyParticleChange.Initialize(track);
  return &fDummyParticleChange;
}


//--------------------------------------------
//         At Rest parameterisation:
//--------------------------------------------
//   AtRestGetPhysiscalInteractionLength:
//--------------------------------------------
G4double 
G4FastSimulationManagerProcess::
AtRestGetPhysicalInteractionLength(const G4Track&        track, 
				   G4ForceCondition* condition)
{
  const G4VPhysicalVolume* currentVolume(0);
  if ( fIsGhostGeometry )  currentVolume = fPathFinder->GetLocatedVolume(fGhostNavigatorIndex);
  else                     currentVolume = track.GetVolume();
  fFastSimulationManager = currentVolume->GetLogicalVolume()->GetFastSimulationManager();
  if( fFastSimulationManager )
    {
      // Ask for trigger:
      fFastSimulationTrigger = fFastSimulationManager->AtRestGetFastSimulationManagerTrigger(track, fGhostNavigator);
      if( fFastSimulationTrigger )
	{
	  // Dirty trick to take control over stepping. Does anyone will ever use that ?
	  *condition = NotForced;
	  return -1.0;
	}
    }
  
  // -- no fast simulation occuring there:
  *condition = NotForced;
  return DBL_MAX;
}

//-----------------------------------------------
//                  AtRestDoIt:
//-----------------------------------------------
G4VParticleChange* G4FastSimulationManagerProcess::AtRestDoIt(const G4Track&, const G4Step&)
{
  return fFastSimulationManager->InvokeAtRestDoIt();
}


void G4FastSimulationManagerProcess::Verbose() const
{
  /*  G4cout << "     >>>>> Trigger Status : ";
  switch(fFastSimulationManager->GetTriggerStatus())
    {
    case NoModel:
      G4cout << "NoModel" << G4endl;
      break;
    case OnBoundaryButLeaving:
      G4cout << "OnBoundaryButLeaving" << G4endl;
      break;
    case OneModelTrigger:
      G4cout << "OneModelTrigger" << G4endl;
      break;
    case NoModelTrigger:
      G4cout << "NoModelTrigger" << G4endl;
      break;
    case Undefined:
      G4cout << "Undefined" << G4endl;
      break;
    default:
      G4cout << " Bizarre..." << G4endl;
      break;
    }*/
}


