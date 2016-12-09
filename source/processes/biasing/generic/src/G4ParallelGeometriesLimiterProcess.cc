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
// $Id:  G4BiasingProcessLimiterForParallelGeometries.cc $
//
//

#include "G4ios.hh"
#include "G4ParallelGeometriesLimiterProcess.hh"
#include "G4BiasingProcessSharedData.hh"
#include "G4ProcessManager.hh"
#include "G4TransportationManager.hh"
#include "G4PathFinder.hh"
#include "G4FieldTrackUpdator.hh"

#include "G4SystemOfUnits.hh"

G4ParallelGeometriesLimiterProcess::G4ParallelGeometriesLimiterProcess(const G4String& processName) :
  G4VProcess(processName, fParallel),
  fParallelWorldSafety( 0.0   ),
  fIsTrackingTime     ( false ),
  fFieldTrack         ( '0'   )
{
  // -- Process Sub Type ? §§
  
  fPathFinder            = G4PathFinder::GetInstance();
  fTransportationManager = G4TransportationManager::GetTransportationManager();
}


// ----------------------------
// -- Add/Remove world volumes:
// ----------------------------
void G4ParallelGeometriesLimiterProcess::AddParallelWorld(const G4String& parallelWorldName)
{
  
  // -- Refuse adding parallel geometry during tracking time:
  if (fIsTrackingTime)
    {
      G4ExceptionDescription ed;
      ed << "G4ParallelGeometriesLimiterProcess `" << GetProcessName()
	 << "': adding a parallel world volume at tracking time is not allowed." << G4endl;
      G4Exception("G4ParallelGeometriesLimiterProcess::AddParallelWorld(const G4String& parallelWorldName)",
		  "BIAS.GEN.21",
		  JustWarning, ed,
		  "Call ignored.");
      return;
    }
  
  else
    
    {
      G4VPhysicalVolume* newWorld = fTransportationManager->IsWorldExisting( parallelWorldName );
      
      // -- Fatal exception if requested world does not exist:
      if (newWorld == 0)
	{
	  G4ExceptionDescription  tellWhatIsWrong;
	  tellWhatIsWrong << "Volume `" <<  parallelWorldName
			  << "' is not a parallel world nor the mass world volume."
			  << G4endl;
	  G4Exception("G4ParallelGeometriesLimiterProcess::SetWorldVolume(const G4String)",
		      "BIAS.GEN.22",
		      FatalException,
		      tellWhatIsWrong);
	}
      
      // -- Protection against adding the mass geometry world as parallel world:
      if ( newWorld ==  fTransportationManager->GetNavigatorForTracking()->GetWorldVolume() )
	{
	  G4ExceptionDescription ed;
	  ed << "G4ParallelGeometriesLimiterProcess `" << GetProcessName()
	     << "': trying to add the world volume for tracking as a parallel world." << G4endl;
	  G4Exception("G4ParallelGeometriesLimiterProcess::AddParallelWorld(const G4String& parallelWorldName)",
		      "BIAS.GEN.23",
		      JustWarning, ed,
		      "Call ignored.");
	  return;
	}
      
      // -- Add parallel world, taking care it is not in the list yet:
      G4bool isNew = true;
      for ( auto knownWorld : fParallelWorlds )
	{
	  if ( knownWorld == newWorld ) isNew = false;
	}
      if ( isNew ) fParallelWorlds.push_back( newWorld );
      else
	{
	  G4ExceptionDescription ed;
	  ed << "G4ParallelGeometriesLimiterProcess `" << GetProcessName()
	     << "': trying to re-add the parallel world volume `" << parallelWorldName << "'." << G4endl;
	  G4Exception("G4ParallelGeometriesLimiterProcess::AddParallelWorld(const G4String& parallelWorldName)",
		      "BIAS.GEN.24",
		      JustWarning, ed,
		      "Call ignored.");
	  return;
	}
    }
  
}


void G4ParallelGeometriesLimiterProcess::RemoveParallelWorld(const G4String& parallelWorldName)
{
  
  // -- Refuse refuse removing parallel geometry during tracking time:
  if (fIsTrackingTime)
    {
      G4ExceptionDescription ed;
      ed << "G4ParallelGeometriesLimiterProcess `" << GetProcessName()
	 << "': removing a parallel world volume at tracking time is not allowed." << G4endl;
      G4Exception("G4ParallelGeometriesLimiterProcess::RemoveParallelWorld(const G4String& parallelWorldName)",
		  "BIAS.GEN.25",
		  JustWarning, ed,
		  "Call ignored.");
      return;
    }
  
  else
    
    {
      G4VPhysicalVolume* newWorld = fTransportationManager->IsWorldExisting( parallelWorldName );
      
      if (newWorld == 0)
	{
	  
	  G4ExceptionDescription ed;
	  ed << "G4ParallelGeometriesLimiterProcess `" << GetProcessName()
	     << "': trying to remove an inexisting parallel world '" << parallelWorldName << "'." << G4endl;
	  G4Exception("G4ParallelGeometriesLimiterProcess::RemoveParallelWorld(const G4String& parallelWorldName)",
		      "BIAS.GEN.26",
		      JustWarning, ed,
		      "Call ignored.");
	  return;
	}
      
      // -- get position of world volume in list:
      size_t iWorld = 0;
      for ( auto knownWorld : fParallelWorlds )
	{
	  if ( knownWorld == newWorld ) break;
	  iWorld++;
	}
      
      if ( iWorld == fParallelWorlds.size() )
	{
	  G4ExceptionDescription ed;
	  ed << "G4ParallelGeometriesLimiterProcess `" << GetProcessName()
	     << "': trying to remove an non-registerered parallel world '" << parallelWorldName << "'." << G4endl;
	  G4Exception("G4ParallelGeometriesLimiterProcess::RemoveParallelWorld(const G4String& parallelWorldName)",
		      "BIAS.GEN.27",
		      JustWarning, ed,
		      "Call ignored.");
	  return;
	}
      
      // -- remove from vector:
      fParallelWorlds.erase( fParallelWorlds.begin() + iWorld );
      
    }
  
  
  
}


// --------------------
//  Start/End tracking:
// --------------------
void G4ParallelGeometriesLimiterProcess::StartTracking(G4Track* track)
{
  fIsTrackingTime = true;
  
  // -- fetch the navigators, their indeces, and activate:
  fParallelWorldNavigators      .clear();
  fParallelWorldNavigatorIndeces.clear();
  fParallelWorldSafeties        .clear();
  fParallelWorldIsLimiting      .clear();
  fParallelWorldWasLimiting     .clear();
  fCurrentVolumes               .clear();
  fPreviousVolumes              .clear();
  for ( auto parallelWorld : fParallelWorlds )
    {
      fParallelWorldNavigators      .push_back( fTransportationManager->     GetNavigator( parallelWorld                   ) );
      fParallelWorldNavigatorIndeces.push_back( fTransportationManager->ActivateNavigator( fParallelWorldNavigators.back() ) );
      fParallelWorldSafeties        .push_back( 0.0 );
      fParallelWorldIsLimiting      .push_back( false );
      fParallelWorldWasLimiting     .push_back( false );
    }
  
  fPathFinder->PrepareNewTrack( track->GetPosition(), track->GetMomentumDirection() );
  // -- §§ does it work at this level, after "PrepareNewTrack" above ?
  for ( auto navigatorIndex : fParallelWorldNavigatorIndeces )
    {
      fPreviousVolumes.push_back( nullptr );
      fCurrentVolumes .push_back( fPathFinder->GetLocatedVolume( navigatorIndex ) );
    }

  // -- will force updating safety:
  fParallelWorldSafety = 0.0;
  for ( size_t i = 0 ; i < fParallelWorldNavigatorIndeces.size() ; i++ ) fParallelWorldSafeties[i] = 0.0;
}


void G4ParallelGeometriesLimiterProcess::EndTracking()
{
  fIsTrackingTime = false;
  for ( auto parallelWorldNavigator : fParallelWorldNavigators )
    fTransportationManager->DeActivateNavigator( parallelWorldNavigator ); 
}


G4double G4ParallelGeometriesLimiterProcess::PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition* condition)
{

  // -- push previous step limitation flags and volumes:
  // -- §§ consider switching pointers insteads of making copies of std::vector's:
  fParallelWorldWasLimiting = fParallelWorldIsLimiting;
  fPreviousVolumes          = fCurrentVolumes;
  
  // -- update volumes:
  size_t i = 0;
  for ( auto navigatorIndex : fParallelWorldNavigatorIndeces ) fCurrentVolumes[i++] = fPathFinder->GetLocatedVolume( navigatorIndex );
  
  *condition = NotForced;
  return DBL_MAX;
}


G4double G4ParallelGeometriesLimiterProcess::AlongStepGetPhysicalInteractionLength(const G4Track&                track,
										   G4double           previousStepSize,
										   G4double         currentMinimumStep, 
										   G4double&            proposedSafety, 
										   G4GPILSelection*          selection)
{
  
  // -- Init:
  // -- Note that the returnedStep must be physically meaningful, even if we return NotCandidateForSelection as condition;
  // -- the reason is that the stepping manager always takes the smallest alongstep among the returned ones (point related
  // -- to geometry step length wrt to true path length).
  *selection            = NotCandidateForSelection;
  G4double returnedStep = DBL_MAX;
  
  // -- G4FieldTrack and ELimited:
  static G4ThreadLocal G4FieldTrack *endTrack_G4MT_TLS_ = 0 ;
  if (!endTrack_G4MT_TLS_) endTrack_G4MT_TLS_ = new  G4FieldTrack ('0') ;
  G4FieldTrack &endTrack = *endTrack_G4MT_TLS_;
  
  static G4ThreadLocal ELimited *eLimited_G4MT_TLS_ = 0 ;
  if (!eLimited_G4MT_TLS_) eLimited_G4MT_TLS_ = new  ELimited  ;
  ELimited &eLimited = *eLimited_G4MT_TLS_;


  // -------------------
  // -- Update safeties:
  // -------------------
  if ( previousStepSize > 0.0 )
    {
      for ( auto& parallelWorldSafety : fParallelWorldSafeties )
	{
	  parallelWorldSafety -= previousStepSize;
	  if ( parallelWorldSafety < 0. ) parallelWorldSafety = 0.0;
	  fParallelWorldSafety =  parallelWorldSafety < fParallelWorldSafety ?  parallelWorldSafety : fParallelWorldSafety ;
	}
    }

  
  // ------------------------------------------
  // Determination of the proposed step length:
  // ------------------------------------------
  if ( ( currentMinimumStep <= fParallelWorldSafety ) && ( currentMinimumStep > 0. ) )
    {
      // -- No chance to limit the step, as proposed move inside safety
      
      returnedStep   = currentMinimumStep;
      proposedSafety = fParallelWorldSafety - currentMinimumStep;
    }
  else
    {
      // -- Proposed move exceeds common safety, need to state
      G4double smallestReturnedStep    = -1.0;
      ELimited eLimitedForSmallestStep = kDoNot;
      for ( size_t i = 0 ; i < fParallelWorldNavigatorIndeces.size() ; i++ )
	{
	  // -- Update safety of geometries having safety smaller than current minimum step
	  if (  currentMinimumStep >= fParallelWorldSafeties[i] )
	    {
	      G4FieldTrackUpdator::Update(&fFieldTrack, &track);
	      G4double tmpReturnedStep = fPathFinder->ComputeStep(fFieldTrack,
								  currentMinimumStep,
								  fParallelWorldNavigatorIndeces[i],
								  track.GetCurrentStepNumber(),
								  fParallelWorldSafeties[i],
								  eLimited,
								  endTrack,
								  track.GetVolume());
	      
	      if ( ( smallestReturnedStep < 0.0 ) || ( tmpReturnedStep <= smallestReturnedStep ) )
		{
		  smallestReturnedStep    = tmpReturnedStep;
		  eLimitedForSmallestStep = eLimited;
		}
	      
	      if (eLimited == kDoNot)
		{
		  // -- Step not limited by this geometry
		  fParallelWorldSafeties[i]   = fParallelWorldNavigators[i]->ComputeSafety(endTrack.GetPosition());
		  fParallelWorldIsLimiting[i] = false;
		}
	      else
		{
		  fParallelWorldIsLimiting[i] = true;
		}
	    }

	  // -- update with smallest safety:
	  fParallelWorldSafety =  fParallelWorldSafeties[i] < fParallelWorldSafety ?  fParallelWorldSafeties[i] : fParallelWorldSafety ;
	}

      // -- no geometry limitation among all geometries, can return currentMinimumStep (or DBL_MAX):
      // -- Beware : the returnedStep must be physically meaningful, even if we say "NotCandidateForSelection" !
      if     (  eLimitedForSmallestStep == kDoNot )
	{
	  returnedStep = currentMinimumStep;
	}
      // -- proposed step length of limiting geometry:
      if     (  eLimitedForSmallestStep == kUnique  ||
	        eLimitedForSmallestStep == kSharedOther )
	{
	  *selection   = CandidateForSelection;
	  returnedStep = smallestReturnedStep;
	}
      else if ( eLimitedForSmallestStep == kSharedTransport)
	{
	  returnedStep = smallestReturnedStep* (1.0 + 1.0e-9);   // -- Expand to disable its selection in Step Manager comparison
	}

      // -- and smallest safety among geometries:
      proposedSafety   = fParallelWorldSafety ;
    }

  // -- returns step length, and proposedSafety
  return returnedStep;
}


G4VParticleChange* G4ParallelGeometriesLimiterProcess::AlongStepDoIt( const G4Track& track,
								      const G4Step&         )
{
  
  fDummyParticleChange.Initialize(track);
  return &fDummyParticleChange;
}


void G4ParallelGeometriesLimiterProcess::SetProcessManager(const G4ProcessManager* mgr)
{
  G4BiasingProcessSharedData *sharedData(nullptr);
  
  // -- initialize sharedData pointer:
  if (  G4BiasingProcessSharedData::fSharedDataMap.Find(mgr) == G4BiasingProcessSharedData::fSharedDataMap.End() )
    {
      sharedData =  new G4BiasingProcessSharedData( mgr );
      G4BiasingProcessSharedData::fSharedDataMap[mgr] = sharedData;
    }
  else sharedData =  G4BiasingProcessSharedData::fSharedDataMap[mgr] ;
  
  // -- add itself to the shared data:
  if ( sharedData->fParallelGeometriesLimiterProcess == nullptr ) sharedData->fParallelGeometriesLimiterProcess = this;
  else
    {
      G4ExceptionDescription ed;
      ed << " Trying to add more than one G4ParallelGeometriesLimiterProcess process to the process manager " << mgr
	 << " (process manager for `" << mgr->GetParticleType()->GetParticleName() << "'). Only one is needed. Call ignored." << G4endl;
      G4Exception("  G4ParallelGeometriesLimiterProcess::SetProcessManager(...)",
		  "BIAS.GEN.29",
		  JustWarning,
		  ed);
    }
}


G4int G4ParallelGeometriesLimiterProcess::GetParallelWorldIndex( const G4VPhysicalVolume* parallelWorld ) const
{
  G4int toReturn = -1;
  G4int iWorld = 0;
  for ( auto world : fParallelWorlds )
    {
      if ( world == parallelWorld )
	{
	  toReturn = iWorld;
	  break;
	}
      iWorld++;
    }
  return toReturn;
}


G4int G4ParallelGeometriesLimiterProcess::GetParallelWorldIndex( G4String parallelWorldName ) const
{
  G4VPhysicalVolume* aWorld = fTransportationManager->IsWorldExisting( parallelWorldName );   // note aWorld might be nullptr
  return GetParallelWorldIndex( aWorld );
}

