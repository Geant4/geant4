//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4SteppingManager.cc,v 1.31 2003-06-16 17:13:19 gunter Exp $
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
  : fUserSteppingAction(NULL), verboseLevel(0)
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
   SetNavigator(G4TransportationManager::GetTransportationManager()
		   ->GetNavigatorForTracking());

   fSelectedAtRestDoItVector
      = new G4SelectedAtRestDoItVector(SizeOfSelectedDoItVector,0);
   fSelectedAlongStepDoItVector
      = new G4SelectedAlongStepDoItVector(SizeOfSelectedDoItVector,0);
   fSelectedPostStepDoItVector
      = new G4SelectedPostStepDoItVector(SizeOfSelectedDoItVector,0);

   SetNavigator(G4TransportationManager::GetTransportationManager()
     ->GetNavigatorForTracking());

   physIntLength = DBL_MAX; 
//   fTouchableHandle = new G4TouchableHistory();
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
   fTrack->SetTouchableHandle(fTrack->GetNextTouchableHandle());

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
     endpointSafety=  std::max( proposedSafety - GeomStepLength, 0.);

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

   PhysicalStep = 0.;
   GeometricalStep = 0.;
   CorrectedStep = 0.;
   PreStepPointIsGeom = false;
   FirstStep = false;
   fStepStatus = fUndefined;

   TempInitVelocity = 0.;
   TempVelocity = 0.;
   sumEnergyChange = 0.;


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


// Set Touchable to track and a private attribute of G4SteppingManager
 

  if ( ! fTrack->GetTouchableHandle() ) {
     G4ThreeVector direction= fTrack->GetMomentumDirection();
     fNavigator->LocateGlobalPointAndSetup( fTrack->GetPosition(), &direction, true, false);
     fTouchableHandle = fNavigator->CreateTouchableHistory();

     fTrack->SetTouchableHandle( fTouchableHandle );
     fTrack->SetNextTouchableHandle( fTouchableHandle );
  }else{
     fTrack->SetNextTouchableHandle( fTrack->GetTouchableHandle() );
     G4VPhysicalVolume* oldTopVolume= fTrack->GetTouchableHandle()->GetVolume();
     G4VPhysicalVolume* newTopVolume=
     fNavigator->LocateGlobalPointAndSetup( fTrack->GetPosition(), 
         fTrack->GetMomentumDirection(),*((G4TouchableHistory*)fTrack->GetTouchableHandle()()) );
     if(newTopVolume != oldTopVolume ){
        fTouchableHandle = fNavigator->CreateTouchableHistory();
        fTrack->SetTouchableHandle( fTouchableHandle );
        fTrack->SetNextTouchableHandle( fTouchableHandle );
     }
  }
// Set vertex information of G4Track at here
   if ( fTrack->GetCurrentStepNumber() == 0 ) {
     fTrack->SetVertexPosition( fTrack->GetPosition() );
     fTrack->SetVertexMomentumDirection( fTrack->GetMomentumDirection() );
     fTrack->SetVertexKineticEnergy( fTrack->GetKineticEnergy() );
     fTrack->SetLogicalVolumeAtVertex( fTrack->GetVolume()->GetLogicalVolume() );
   }
// Initial set up for attributes of 'G4SteppingManager'
   fCurrentVolume = fTouchableHandle->GetVolume();

// If track is already outside the world boundary, kill it
   if( fCurrentVolume==0 ){
       // If the track is a primary, stop processing
       if(fTrack->GetParentID()==0)
       {
        G4cerr << "Primary particle starting at "
               << fTrack->GetPosition()
               << " is outside of the world volume." << G4endl;
        G4Exception("G4SteppingManager::Primary vertex outside of the world");
       }

       fTrack->SetTrackStatus( fStopAndKill );
       G4cerr << "G4SteppingManager::SetInitialStep(): warning: "
              << "initial track position is outside world! "
              << fTrack->GetPosition() << G4endl;
   }
   else {
// Initial set up for attribues of 'Step'
       fStep->InitializeStep( fTrack );
   }
#ifdef G4VERBOSE
                         // !!!!! Verbose
           if(verboseLevel>0) fVerbose->TrackingStarted();
#endif
}

