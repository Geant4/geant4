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
// $Id: G4ErrorPropagator.cc 91809 2015-08-06 13:00:09Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorPropagator.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorFreeTrajState.hh"
#include "G4ErrorSurfaceTrajState.hh"
#include "G4ErrorGeomVolumeTarget.hh"
#include "G4ErrorSurfaceTarget.hh"

#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4ParticleTable.hh"
#include "G4StateManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4UnitsTable.hh"

#include <vector>


//---------------------------------------------------------------------------
G4ErrorPropagator::G4ErrorPropagator()
  : theStepLength(0.), theInitialTrajState(0), theStepN(0), theG4Track(0)
{
  verbose =  G4ErrorPropagatorData::verbose();
#ifdef G4EVERBOSE
   if(verbose >= 5) { G4cout << "G4ErrorPropagator " << this << G4endl; }
#endif

  fpSteppingManager = G4EventManager::GetEventManager()
                    ->GetTrackingManager()->GetSteppingManager();
  thePropIsInitialized = false;
}


//-----------------------------------------------------------------------
G4int G4ErrorPropagator::Propagate( G4ErrorTrajState* currentTS,
                              const G4ErrorTarget* target, G4ErrorMode mode )
{
  // to start ierror is set to 1 (= OK)
  //
  G4int ierr = 1;

  G4ErrorPropagatorData* g4edata =
    G4ErrorPropagatorData::GetErrorPropagatorData();

  //----- Do not propagate zero or too low energy particles
  //
  if( currentTS->GetMomentum().mag() < 1.E-9*MeV )
  {
    std::ostringstream message;
    message << "Energy too low to be propagated: "
            << G4BestUnit(currentTS->GetMomentum().mag(),"Energy");
    G4Exception("G4ErrorPropagator::Propagate()", "GEANT4e-Notification",
                JustWarning, message);
    return -3; 
  }

  g4edata->SetMode( mode );

#ifdef G4EVERBOSE
  if( verbose >= 1 )
  {
     G4cout << " =====> starting GEANT4E tracking for "
            << currentTS->GetParticleType()
            << "  Forwards= " << g4edata->GetMode() << G4endl;
  }
  if(verbose >= 1 )
  {
     G4cout << G4endl << "@@@@@@@@@@@@@@@@@@@@@@@@@ NEW STEP " << G4endl;
  }

  if( verbose >= 3 )
  {
    G4cout << " G4ErrorPropagator::Propagate initialTS ";
    G4cout << *currentTS << G4endl;
    target->Dump(G4String(" to target "));
  }
#endif

  g4edata->SetTarget( target );

  //----- Create a track
  //
  if( theG4Track != 0 ) { delete theG4Track; }
  theG4Track = InitG4Track( *currentTS );

  //----- Create a G4ErrorFreeTrajState
  //
  G4ErrorFreeTrajState* currentTS_FREE = InitFreeTrajState( currentTS );

  //----- Track the particle
  //
  ierr = MakeSteps( currentTS_FREE );

  //------ Tracking ended, check if target has been reached
  //       if target not found
  //
  if( g4edata->GetState() != G4ErrorState_StoppedAtTarget )
  {
    if( theG4Track->GetKineticEnergy() > 0. )
    {
      ierr = -ierr - 10;
    }
    else
    {
      ierr = -ierr - 20;
    }
    *currentTS = *currentTS_FREE;
    if(verbose >= 0 )
    {
      std::ostringstream message;
      message << "Particle does not reach target: " << *currentTS;
      G4Exception("G4ErrorPropagator::Propagate()", "GEANT4e-Notification",
                  JustWarning, message);
    }
  }
  else
  {
    GetFinalTrajState( currentTS, currentTS_FREE, target );
  }

#ifdef G4EVERBOSE
  if( verbose >= 1 )
  {
    G4cout << " G4ErrorPropagator: propagation ended " << G4endl;
  }
  if( verbose >= 2 )
  {
    G4cout << " Current TrajState " << currentTS << G4endl;
  }
#endif
 
  // Inform end of tracking to physics processes
  //
  theG4Track->GetDefinition()->GetProcessManager()->EndTracking();

  InvokePostUserTrackingAction( theG4Track );

  // delete currentTS_FREE;

  return ierr;
}


//-----------------------------------------------------------------------
G4int G4ErrorPropagator::PropagateOneStep( G4ErrorTrajState* currentTS )
{
  G4ErrorPropagatorData* g4edata =
    G4ErrorPropagatorData::GetErrorPropagatorData();

  if ( (g4edata->GetState() == G4ErrorState_PreInit)
    || (G4StateManager::GetStateManager()->GetCurrentState()
       != G4State_GeomClosed) )
  {
    std::ostringstream message;
    message << "Called before initialization is done for this track!";
    G4Exception("G4ErrorPropagator::PropagateOneStep()",
                "InvalidCall", FatalException, message,
                "Please call G4ErrorPropagatorManager::InitGeant4e().");
  }

  // to start ierror is set to 0 (= OK)
  //
  G4int ierr = 0;

  //--- Do not propagate zero or too low energy particles
  //
  if( currentTS->GetMomentum().mag() < 1.E-9*MeV )
  {
    std::ostringstream message;
    message << "Energy too low to be propagated: "
            << G4BestUnit(currentTS->GetMomentum().mag(),"Energy");
    G4Exception("G4ErrorPropagator::PropagateOneStep()",
                "GEANT4e-Notification", JustWarning, message);
    return -3;   
  }

#ifdef G4EVERBOSE
  if( verbose >= 1 )
  {
    G4cout << " =====> starting GEANT4E tracking for "
           << currentTS->GetParticleType()
           << "  Forwards= " << g4edata->GetMode() << G4endl;
  }
  if( verbose >= 3 )
  {
    G4cout << " G4ErrorPropagator::Propagate initialTS ";
    G4cout << *currentTS << G4endl;
  }
#endif

  //----- If it is the first step, create a track
  //
  if( theStepN == 0 )
  {
    if( theG4Track != 0 ) { delete theG4Track; }
    theG4Track = InitG4Track( *currentTS );
  }
    // set to 0 by the initialization in G4ErrorPropagatorManager
  theStepN++;

  //----- Create a G4ErrorFreeTrajState
  //
  G4ErrorFreeTrajState* currentTS_FREE = InitFreeTrajState( currentTS );

  //----- Track the particle one step
  //
  ierr = MakeOneStep( currentTS_FREE );

  //----- Get the state on target
  //
  GetFinalTrajState( currentTS, currentTS_FREE, g4edata->GetTarget() );

  return ierr;
}


//---------------------------------------------------------------------------
G4Track* G4ErrorPropagator::InitG4Track( G4ErrorTrajState& initialTS )
{
  if( verbose >= 5 ) { G4cout << "InitG4Track " << G4endl; }

  //----- Create Particle
  //
  const G4String partType = initialTS.GetParticleType();
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(partType); 
  if( particle == 0)
  {
    std::ostringstream message;
    message << "Particle type not defined: " << partType;
    G4Exception("G4ErrorPropagator::InitG4Track()", "InvalidSetup",
                FatalException, message);
  }
 
  G4DynamicParticle* DP = 
    new G4DynamicParticle(particle,initialTS.GetMomentum());

  DP->SetPolarization(0.,0.,0.);

  // Set Charge
  //
  if( particle->GetPDGCharge() < 0 )
  {
    DP->SetCharge(-eplus);
  }
  else
  {
    DP->SetCharge(eplus);
  }

  //----- Create Track 
  //
  theG4Track = new G4Track(DP, 0., initialTS.GetPosition() );
  theG4Track->SetParentID(0);

#ifdef G4EVERBOSE
  if(verbose >= 3)
  {
    G4cout << " G4ErrorPropagator new track of energy: "
           << theG4Track->GetKineticEnergy() << G4endl;
  }
#endif
  
  //---- Reproduce G4TrackingManager::ProcessOneTrack initialization
  InvokePreUserTrackingAction( theG4Track );  

  if( fpSteppingManager == 0 )
  {
    G4Exception("G4ErrorPropagator::InitG4Track()", "InvalidSetup",
                FatalException, "G4SteppingManager not initialized yet!");
  }
  else
  {
    fpSteppingManager->SetInitialStep(theG4Track);
  }

  // Give SteppingManger the maximum number of processes
  //
  fpSteppingManager->GetProcessNumber();

  // Give track the pointer to the Step
  //
  theG4Track->SetStep(fpSteppingManager->GetStep());

  // Inform beginning of tracking to physics processes
  //
  theG4Track->GetDefinition()->GetProcessManager()->StartTracking(theG4Track);

  initialTS.SetG4Track( theG4Track );

  return theG4Track;
}


//-----------------------------------------------------------------------
G4int G4ErrorPropagator::MakeSteps( G4ErrorFreeTrajState* currentTS_FREE )
{
  G4int ierr = 0;

  //----- Track the particle Step-by-Step while it is alive
  //
  theStepLength = 0.;
  
  while( (theG4Track->GetTrackStatus() == fAlive) ||
         (theG4Track->GetTrackStatus() == fStopButAlive) )
  {
    ierr = MakeOneStep( currentTS_FREE );
    if( ierr != 0 ) { break; }

    //----- Check if last step for error propagation
    //
    if( CheckIfLastStep( theG4Track ) )
    {
      break;
    }
  }  // Loop checking, 06.08.2015, G.Cosmo
  return ierr;
}


//-----------------------------------------------------------------------
G4int G4ErrorPropagator::MakeOneStep( G4ErrorFreeTrajState* currentTS_FREE )
{
  G4ErrorPropagatorData* g4edata =
    G4ErrorPropagatorData::GetErrorPropagatorData();
  G4int ierr = 0;

  //---------- Track one step
#ifdef G4EVERBOSE
  if(verbose >= 2 )
  {
    G4cout << G4endl
           << "@@@@@@@@@@@@@@@@@@@@@@@@@ NEW STEP " << G4endl;
  }
#endif
  
  theG4Track->IncrementCurrentStepNumber();

  fpSteppingManager->Stepping();
  
  //---------- Check if Target has been reached (and then set G4ErrorState)

  // G4ErrorPropagationNavigator limits the step if target is closer than
  // boundary (but the winner process is always "Transportation": then
  // error propagator will stop the track

  if( theG4Track->GetStep()->GetPostStepPoint()
      ->GetProcessDefinedStep()->GetProcessName() == "Transportation" )
  {
    if( g4edata->GetState()
       == G4ErrorState(G4ErrorState_TargetCloserThanBoundary) )
    {  // target or step length reached
      
#ifdef G4EVERBOSE
      if(verbose >= 5 )
      {
        G4cout << " transportation determined by geant4e " << G4endl;
      }
#endif      
      g4edata->SetState( G4ErrorState_StoppedAtTarget );
    }
    else if( g4edata->GetTarget()->GetType() == G4ErrorTarget_GeomVolume )
    {
      G4ErrorGeomVolumeTarget* target =
        (G4ErrorGeomVolumeTarget*)(g4edata->GetTarget());
      if( static_cast<G4ErrorGeomVolumeTarget*>( target )
          ->TargetReached( theG4Track->GetStep() ) )
      {
        g4edata->SetState( G4ErrorState_StoppedAtTarget ); 
      } 
    }
  }
  else if( theG4Track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()
           ->GetProcessName() == "G4ErrorTrackLengthTarget" )
  {
    g4edata->SetState( G4ErrorState_StoppedAtTarget );
  }

  //---------- Propagate error  

#ifdef G4EVERBOSE
  if(verbose >= 2 )
  {
    G4cout << " propagating error " << *currentTS_FREE << G4endl;
  }
#endif
  const G4Track* cTrack = const_cast<G4Track*>(theG4Track);
  ierr = currentTS_FREE->PropagateError( cTrack );
  
#ifdef G4EVERBOSE
  if(verbose >= 3 )
  {
    G4cout << " PropagateError returns " << ierr << G4endl;
  }
#endif

  currentTS_FREE->Update( cTrack );
  
  theStepLength += theG4Track->GetStepLength();
   
  if(ierr != 0 )
  {
    std::ostringstream message;
    message << "Error returned: " << ierr;
    G4Exception("G4ErrorPropagator::MakeOneStep()",
                "GEANT4e-Notification", JustWarning, message,
                "Geant4 tracking will be stopped !");
  }

  return ierr; 
}


//-----------------------------------------------------------------------
G4ErrorFreeTrajState*
G4ErrorPropagator::InitFreeTrajState( G4ErrorTrajState* currentTS )
{
  G4ErrorFreeTrajState* currentTS_FREE = 0;

  //----- Transform the TrajState to Free coordinates if it is OnSurface
  //
  if( currentTS->GetTSType() == G4eTS_OS )
  {
    G4ErrorSurfaceTrajState* tssd =
      static_cast<G4ErrorSurfaceTrajState*>(currentTS);
    currentTS_FREE = new G4ErrorFreeTrajState( *tssd );
  }
  else if( currentTS->GetTSType() == G4eTS_FREE )
  {
    currentTS_FREE = static_cast<G4ErrorFreeTrajState*>(currentTS);
  }
  else
  {
    std::ostringstream message;
    message << "Wrong trajectory state: " << currentTS->GetTSType();
    G4Exception("G4ErrorPropagator::InitFreeTrajState()", "InvalidState",
                FatalException, message);
  }
  return currentTS_FREE;
}


//-----------------------------------------------------------------------
void G4ErrorPropagator::GetFinalTrajState( G4ErrorTrajState* currentTS,
                                           G4ErrorFreeTrajState* currentTS_FREE,
                                           const G4ErrorTarget* target ) 
{
  G4ErrorPropagatorData* g4edata =
    G4ErrorPropagatorData::GetErrorPropagatorData();

#ifdef G4EVERBOSE
  if(verbose >= 1 )
  {
    G4cout << " G4ErrorPropagator::Propagate: final state "
           << G4int(g4edata->GetState()) << " TSType "
           << currentTS->GetTSType() << G4endl;
  }
#endif

  if( (currentTS->GetTSType() == G4eTS_FREE) || 
      (g4edata->GetState() != G4ErrorState_StoppedAtTarget) )
  {
    currentTS = currentTS_FREE;
  }
  else if( currentTS->GetTSType() == G4eTS_OS )
  {
    if( target->GetType() == G4ErrorTarget_TrkL )
    {
      G4Exception("G4ErrorPropagator:GetFinalTrajState()",
                  "InvalidSetup", FatalException,
                  "Using a G4ErrorSurfaceTrajState with wrong target");
    }
    const G4ErrorTanPlaneTarget* targetWTP =
      static_cast<const G4ErrorTanPlaneTarget*>(target);
    *currentTS = G4ErrorSurfaceTrajState(
                 *(static_cast<G4ErrorFreeTrajState*>(currentTS_FREE)),
                 targetWTP->GetTangentPlane( currentTS_FREE->GetPosition() ) );
#ifdef G4EVERBOSE
    if(verbose >= 1 )
    {
      G4cout << currentTS << " returning tssd " << *currentTS << G4endl;
    }
#endif
    delete currentTS_FREE;
  }
}


//-------------------------------------------------------------------------
G4bool G4ErrorPropagator::CheckIfLastStep( G4Track* aTrack )
{
  G4bool exception = 0;
  G4bool lastG4eStep = false;
  G4ErrorPropagatorData* g4edata =
    G4ErrorPropagatorData::GetErrorPropagatorData();

#ifdef G4EVERBOSE
  if( verbose >= 4 )
  {
    G4cout << " G4ErrorPropagator::CheckIfLastStep G4ErrorState= "
           << G4int(g4edata->GetState()) << G4endl;
  }
#endif
  
  //----- Check if this is the last step: track has reached the target
  //      or the end of world
  //
  if(g4edata->GetState() == G4ErrorState(G4ErrorState_StoppedAtTarget) )
  {
    lastG4eStep = true;    
#ifdef G4EVERBOSE
    if(verbose >= 5 )
    {
      G4cout << " G4ErrorPropagator::CheckIfLastStep " << lastG4eStep
             << " " << G4int(g4edata->GetState()) << G4endl;
    }
#endif
  }
  else if( aTrack->GetNextVolume() == 0 )
  {
    //----- If particle is out of world, without finding the G4ErrorTarget,
    //      give a n error/warning
    //
    lastG4eStep = true;
    if( exception )
    {
      std::ostringstream message;
      message << "Track extrapolated until end of World" << G4endl
              << "without finding the defined target!";
      G4Exception("G4ErrorPropagator::CheckIfLastStep()",
                  "InvalidSetup", FatalException, message);
    }
    else
    {
      if( verbose >= 1 )
      {
        std::ostringstream message;
        message << "Track extrapolated until end of World" << G4endl
                << "without finding the defined target.";
        G4Exception("G4ErrorPropagator::CheckIfLastStep()",
                    "GEANT4e-Notification", JustWarning, message);
      }
    }
  }  //----- not last step from G4e, but track is stopped (energy exhausted)
  else if( aTrack->GetTrackStatus() == fStopAndKill )
  { 
    if( exception )
    {
      std::ostringstream message;
      message << "Track extrapolated until energy is exhausted" << G4endl
              << "without finding the defined target!";
      G4Exception("G4ErrorPropagator::CheckIfLastStep()",
                  "InvalidSetup", FatalException, message);
    }
    else
    {
      if( verbose >= 1 )
      {
        std::ostringstream message;
        message << "Track extrapolated until energy is exhausted" << G4endl
                << "without finding the defined target.";
        G4Exception("G4ErrorPropagator::CheckIfLastStep()",
                    "GEANT4e-Notification", JustWarning, message);
      }
      lastG4eStep = 1;
    }
  }

#ifdef G4EVERBOSE
  if( verbose >= 5 )
  {
    G4cout << " return CheckIfLastStep " << lastG4eStep << G4endl;
  }
#endif

  return  lastG4eStep;
}


//---------------------------------------------------------------------------
void G4ErrorPropagator::InvokePreUserTrackingAction( G4Track* fpTrack )
{
  const G4UserTrackingAction* fpUserTrackingAction =
    G4EventManager::GetEventManager()->GetUserTrackingAction();
  if( fpUserTrackingAction != 0 )
  {
    const_cast<G4UserTrackingAction*>(fpUserTrackingAction)
      ->PreUserTrackingAction((fpTrack) );
  }
}


//---------------------------------------------------------------------------
void G4ErrorPropagator::InvokePostUserTrackingAction( G4Track* fpTrack )
{
  const G4UserTrackingAction* fpUserTrackingAction =
    G4EventManager::GetEventManager()->GetUserTrackingAction();
  if( fpUserTrackingAction != 0 )
  {
    const_cast<G4UserTrackingAction*>(fpUserTrackingAction)
      ->PostUserTrackingAction((fpTrack) );
  }
}
