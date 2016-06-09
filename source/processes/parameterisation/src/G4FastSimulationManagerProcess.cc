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
// $Id: G4FastSimulationManagerProcess.cc,v 1.15 2007/05/11 13:50:20 mverderi Exp $
// GEANT4 tag $Name: geant4-09-00 $
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

#define  PARANOIA

G4FastSimulationManagerProcess::
G4FastSimulationManagerProcess(const G4String& processName,
							G4ProcessType       theType) : 
  G4VProcess(processName,theType),
  fWorldVolume(0),
  fIsTrackingTime(false),
  fNavigator(0),
  fNavigatorIndex(-1),
  fFastSimulationManager(0),
  fFastSimulationTrigger(false)
{
  SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume()->GetName());
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
  fWorldVolume(0),
  fIsTrackingTime(false),
  fNavigator(0),
  fNavigatorIndex(-1),
  fFastSimulationManager(0),
  fFastSimulationTrigger(false)
{
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
  fWorldVolume(0),
  fIsTrackingTime(false),
  fNavigator(0),
  fNavigatorIndex(-1),
  fFastSimulationManager(0),
  fFastSimulationTrigger(false)
{
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
void
G4FastSimulationManagerProcess::
SetWorldVolume(G4String newWorldName)
{
  if (fIsTrackingTime)
    G4cout << "!!! G4FastSimulationManagerProcess `" << GetProcessName() 
	   << "': changing world volume at tracking time is not allowed for now. Call ignored !!!" << G4endl;
  else
    {
      G4VPhysicalVolume* newWorld = G4TransportationManager::GetTransportationManager()->IsWorldExisting(newWorldName);
      G4String tellWhatIsWrong;
      tellWhatIsWrong = "Volume newWorldName = `"; tellWhatIsWrong +=  newWorldName; tellWhatIsWrong += "' is not a parallel world nor the mass world volume.";
      if (newWorld == 0) G4Exception("G4FastSimulationManagerProcess::SetWorldVolume(const G4String&, G4bool verbose)",
				     "InvalidWorld",
				     FatalException,
				     tellWhatIsWrong);
      if (verboseLevel>0)
	if (fWorldVolume) G4cout << "G4FastSimulationManagerProcess `" << GetProcessName()
				 << "': changing world volume from '"  << fWorldVolume->GetName() 
				 << "' to `" << newWorld << "'." << G4endl;
	else              G4cout << "G4FastSimulationManagerProcess `" << GetProcessName()
				 << "': setting world volume from to `"<< newWorld->GetName() << "'." << G4endl;
      fWorldVolume = newWorld;
    }
  
}


void
G4FastSimulationManagerProcess::
SetWorldVolume(G4VPhysicalVolume* newWorld)
{
  if (newWorld) SetWorldVolume(newWorld->GetName());
  else G4cout << "!!! G4FastSimulationManagerProcess::SetWorldVolume(const G4VPhysicalVolume* newWorld) : null pointer passed. !!! Continuing at your own risks..." << G4endl;
}


// --------------------
//  Start/End tracking:
// --------------------
void
G4FastSimulationManagerProcess::
StartTracking(G4Track*)
{
  fIsTrackingTime = true;
  fIsFirstStep    = true;
  
  // -- fetch the navigator (and its index) and activate it:
  G4TransportationManager* transportationManager = G4TransportationManager::GetTransportationManager();
  fNavigator = transportationManager->GetNavigator(fWorldVolume);
  if (fNavigator != transportationManager->GetNavigatorForTracking())
    fNavigatorIndex = transportationManager->ActivateNavigator(fNavigator);
  else
    fNavigatorIndex = 0;
}


void
G4FastSimulationManagerProcess::
EndTracking()
{
  fIsTrackingTime = false;
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
#ifdef PARANOIA
  if ( fNavigator->GetWorldVolume() != fWorldVolume ) G4Exception("!!! ??? INCONSISTENT NAVIGATORS/WORLD VOLUMES ??? !!!");
#endif
  // -- Get current volume, and check for presence of fast simulation manager.
  // -- For the case of the navigator for tracking (fNavigatorIndex == 0)
  // -- we use the track volume. This allows the code to be valid for both
  // -- cases where the PathFinder is used (G4CoupledTranportation) or not
  // -- (G4Transportation).
  const G4VPhysicalVolume* currentVolume(0);
  if (fNavigatorIndex == 0) currentVolume = track.GetVolume();
  else                      currentVolume = G4PathFinder::GetInstance()->GetLocatedVolume(fNavigatorIndex);
  if ( currentVolume )
    {
      fFastSimulationManager = currentVolume->GetLogicalVolume()->GetFastSimulationManager();
      if( fFastSimulationManager )
	{
	  // Ask for trigger:
	  fFastSimulationTrigger = fFastSimulationManager->PostStepGetFastSimulationManagerTrigger(track, fNavigator);
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


//--------------------------------------------
//         At Rest parameterisation:
//--------------------------------------------
//   AtRestGetPhysiscalInteractionLength:
//--------------------------------------------
G4double 
G4FastSimulationManagerProcess::
AtRestGetPhysicalInteractionLength(const G4Track&    track, 
				   G4ForceCondition* condition)
{
  const G4VPhysicalVolume* currentVolume(0);
  if (fNavigatorIndex == 0) currentVolume = track.GetVolume();
  else                      currentVolume = G4PathFinder::GetInstance()->GetLocatedVolume(fNavigatorIndex);
  fFastSimulationManager = currentVolume->GetLogicalVolume()->GetFastSimulationManager();
  if( fFastSimulationManager )
    {
      // Ask for trigger:
      fFastSimulationTrigger = fFastSimulationManager->AtRestGetFastSimulationManagerTrigger(track, fNavigator);
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


