// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SteppingManager.cc,v 1.7 1999-10-22 03:21:47 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4SteppingManager.cc
//
// Description:
//   This class represents the manager who steers to move the give
//   particle from the TrackingManger by one Step.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#include "G4SteppingManager.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4SteppingControl.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4VSensitiveDetector.hh"    // Include from 'hits/digi'

//////////////////////////////////////
G4SteppingManager::G4SteppingManager()
//////////////////////////////////////
  : verboseLevel(0), fUserSteppingAction(NULL)
{

// Construct simple 'has-a' related objects
   fSecondary = new G4TrackVector();
   fStep = new G4Step();
   fPreStepPoint  = fStep->GetPreStepPoint();
   fPostStepPoint = fStep->GetPostStepPoint();
#ifdef G4VERBOSE
   if(G4VSteppingVerbose::GetInstance()==0) {
     fVerbose =  new G4SteppingVerbose();
     G4VSteppingVerbose::SetInstance(fVerbose);
     fVerbose -> SetManager(this);
     KillVerbose = true;
   }
   else { 
      fVerbose = G4VSteppingVerbose::GetInstance();
      fVerbose -> SetManager(this);
	  KillVerbose = false;
   }
#endif
   fSelectedAtRestDoItVector 
      = new G4SelectedAtRestDoItVector(SIZEofSelectedDoIt);
   fSelectedAlongStepDoItVector 
      = new G4SelectedAlongStepDoItVector(SIZEofSelectedDoIt);
   fSelectedPostStepDoItVector 
      = new G4SelectedPostStepDoItVector(SIZEofSelectedDoIt);
   SetNavigator(G4TransportationManager::GetTransportationManager()
     ->GetNavigatorForTracking());

   fTouchable1 = new G4TouchableHistory();
   fTouchable2 = new G4TouchableHistory();
}

///////////////////////////////////////
G4SteppingManager::~G4SteppingManager()
///////////////////////////////////////
{

// Destruct simple 'has-a' objects
   delete fSecondary;
   delete fStep;
   delete fSelectedAtRestDoItVector;
   delete fSelectedAlongStepDoItVector;
   delete fSelectedPostStepDoItVector;
   delete fTouchable1;
   delete fTouchable2;
   if (fUserSteppingAction) delete fUserSteppingAction;
#ifdef G4VERBOSE
   if(KillVerbose) delete fVerbose;
#endif
}


//////////////////////////////////////////
G4StepStatus G4SteppingManager::Stepping()
//////////////////////////////////////////
{

//--------
// Prelude
//--------
#ifdef G4VERBOSE
            // !!!!! Verbose
             if(verboseLevel>0) fVerbose->NewStep();
#endif 

// Store last PostStepPoint to PreStepPoint, and swap current and nex
// volume information of G4Track. Reset total energy deposit in one Step. 
   fStep->CopyPostToPreStepPoint();
   fStep->ResetTotalEnergyDeposit();

// Switch next touchable in track to current one
   fTrack->SetTouchable(fTrack->GetNextTouchable());

// Free the touchable which is not used
   SetAnotherTouchableFree(fPostStepPoint->GetTouchable());

// Reset the secondary particles
   fN2ndariesAtRestDoIt = 0;
   fN2ndariesAlongStepDoIt = 0;
   fN2ndariesPostStepDoIt = 0;

//JA Set the volume before it is used (in DefineStepLength() for User Limit) 
   fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();

//-----------------
// AtRest Processes
//-----------------

   if( fTrack->GetTrackStatus() == fStopButAlive ){
     if( MAXofAtRestLoops>0 ){
        InvokeAtRestDoItProcs();
        fStepStatus = fAtRestDoItProc;
        fStep->GetPostStepPoint()->SetStepStatus( fStepStatus );
       
#ifdef G4VERBOSE
            // !!!!! Verbose
             if(verboseLevel>0) fVerbose->AtRestDoItInvoked();
#endif 

     }
     // Make sure the track is killed
     fTrack->SetTrackStatus( fStopAndKill );
   }

//---------------------------------
// AlongStep and PostStep Processes
//---------------------------------


   else{
     // Find minimum Step length demanded by active disc./cont. processes
     DefinePhysicalStepLength();

     // Store the Step length (geometrical length) to G4Step and G4Track
     fStep->SetStepLength( PhysicalStep );
     fTrack->SetStepLength( PhysicalStep );
     G4double GeomStepLength = PhysicalStep;

     // Store StepStatus to PostStepPoint
     fStep->GetPostStepPoint()->SetStepStatus( fStepStatus );

     // Invoke AlongStepDoIt 
     InvokeAlongStepDoItProcs();

     // Update track by taking into account all changes by AlongStepDoIt
     fStep->UpdateTrack();

     // Update safety after invocation of all AlongStepDoIts
     endpointSafOrigin= fPostStepPoint->GetPosition();
     endpointSafety=  max( proposedSafety - GeomStepLength, 0.);

     fStep->GetPostStepPoint()->SetSafety( endpointSafety );

#ifdef G4VERBOSE
                         // !!!!! Verbose
           if(verboseLevel>0) fVerbose->AlongStepDoItAllDone();
#endif

     // Invoke PostStepDoIt
     InvokePostStepDoItProcs();

#ifdef G4VERBOSE
                 // !!!!! Verbose
     if(verboseLevel>0) fVerbose->PostStepDoItAllDone();
#endif
   }

//-------
// Finale
//-------

// Update 'TrackLength' and remeber the Step length of the current Step
   fTrack->AddTrackLength(fStep->GetStepLength());
   fPreviousStepSize = fStep->GetStepLength();
#ifdef G4VERBOSE
                         // !!!!! Verbose
           if(verboseLevel>0) fVerbose->StepInfo();
#endif
// Send G4Step information to Hit/Dig if the volume is sensitive
   fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();
   StepControlFlag =  fStep->GetControlFlag();
   if( fCurrentVolume != 0 && StepControlFlag != AvoidHitInvocation) {
      fSensitive = fCurrentVolume->GetLogicalVolume()->
                                   GetSensitiveDetector();
      if( fSensitive != 0 ) {
        fSensitive->Hit(fStep);
      }
   }

// User intervention process.
   fStep->SetTrack(fTrack);
   if( fUserSteppingAction != NULL ) {
      fUserSteppingAction->UserSteppingAction(fStep);
   }

// Stepping process finish. Return the value of the StepStatus.
   return fStepStatus;

}

///////////////////////////////////////////////////////////
void G4SteppingManager::SetInitialStep(G4Track* valueTrack)
///////////////////////////////////////////////////////////
{

// Set up several local variables.
   PreStepPointIsGeom = false;
   FirstStep = true;
   fParticleChange = NULL;
   fPreviousStepSize = 0.;
   fStepStatus = fUndefined;

   fTrack = valueTrack;
   Mass = fTrack->GetDynamicParticle()->GetMass();

   fIsTouchable1Free = TRUE;
   fIsTouchable2Free = TRUE;

// If the primary track has 'Suspend' or 'PostponeToNextEvent' state,
// set the track state to 'Alive'.
   if( (fTrack->GetTrackStatus()==fSuspend) ||
       (fTrack->GetTrackStatus()==fPostponeToNextEvent) ){ 
      fTrack->SetTrackStatus(fAlive);
   }

// If the primary track has 'zero' kinetic energy, set the track
// state to 'StopButAlive'.
   if(fTrack->GetKineticEnergy() <= 0.0){
     fTrack->SetTrackStatus( fStopButAlive );
   }

// Initialize VTouchable using the point in the global coordinate
// system. 
   G4VTouchable* pTouchableFree = GetFreeTouchable();
   fNavigator->LocateGlobalPointAndUpdateTouchable(
                     fTrack->GetPosition(),
                     pTouchableFree,
                     true );

// Set Touchable to track and a private attribute of G4SteppingManager
   fTrack->SetTouchable( pTouchableFree );
   fTrack->SetNextTouchable( pTouchableFree );

// Set vertex information of G4Track at here
   fTrack->SetVertexPosition( fTrack->GetPosition() );
   fTrack->SetVertexMomentumDirection( fTrack->GetMomentumDirection() );
   fTrack->SetVertexKineticEnergy( fTrack->GetKineticEnergy() );

// Initial set up for attributes of 'G4SteppingManager'
   fCurrentVolume = pTouchableFree->GetVolume();

// Initial set up for attribues of 'Step'
   fStep->InitializeStep( fTrack );
#ifdef G4VERBOSE
                         // !!!!! Verbose
           if(verboseLevel>0) fVerbose->TrackingStarted();
#endif
// If track is already outside the world boundary, kill it
   if( fTrack->GetVolume()==0 ){
       fTrack->SetTrackStatus( fStopAndKill );
   }
}

