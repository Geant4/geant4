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
// $Id: G4SteppingManager.cc 95198 2016-01-29 08:28:27Z gcosmo $
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
#include "G4GeometryTolerance.hh"

//////////////////////////////////////
G4SteppingManager::G4SteppingManager()
//////////////////////////////////////
  : fUserSteppingAction(0), verboseLevel(0)
{

// Construct simple 'has-a' related objects
   fStep = new G4Step();
   fSecondary = fStep->NewSecondaryVector();
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
   kCarTolerance = 0.5*G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

///////////////////////////////////////
G4SteppingManager::~G4SteppingManager()
///////////////////////////////////////
{
   fTouchableHandle = 0;
// Destruct simple 'has-a' objects
   fStep->DeleteSecondaryVector();
///////////////////////////   delete fSecondary;
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
	     else 
             if(verboseLevel==-1) { 
		 G4VSteppingVerbose::SetSilent(1);
	     }
	     else
 		 G4VSteppingVerbose::SetSilent(0);
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

// Reset the step's auxiliary points vector pointer
   fStep->SetPointerToVectorOfAuxiliaryPoints(0);

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
//     endpointSafety=  std::max( proposedSafety - GeomStepLength, 0.);
     endpointSafety=  std::max( proposedSafety - GeomStepLength, kCarTolerance);

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
   fStep->SetTrack(fTrack);
#ifdef G4VERBOSE
                         // !!!!! Verbose

           if(verboseLevel>0) fVerbose->StepInfo();
#endif
// Send G4Step information to Hit/Dig if the volume is sensitive
   fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();
   StepControlFlag =  fStep->GetControlFlag();
   if( fCurrentVolume != 0 && StepControlFlag != AvoidHitInvocation) {
      fSensitive = fStep->GetPreStepPoint()->
                                   GetSensitiveDetector();
      if( fSensitive != 0 ) {
        fSensitive->Hit(fStep);
      }
   }

// User intervention process.
   if( fUserSteppingAction != 0 ) {
      fUserSteppingAction->UserSteppingAction(fStep);
   }
   G4UserSteppingAction* regionalAction
    = fStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetRegion()
      ->GetRegionalSteppingAction();
   if( regionalAction ) regionalAction->UserSteppingAction(fStep);

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
   fParticleChange = 0;
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
     fNavigator->LocateGlobalPointAndSetup( fTrack->GetPosition(),
                                            &direction, false, false );
     fTouchableHandle = fNavigator->CreateTouchableHistory();

     fTrack->SetTouchableHandle( fTouchableHandle );
     fTrack->SetNextTouchableHandle( fTouchableHandle );
  }else{
     fTrack->SetNextTouchableHandle( fTouchableHandle = fTrack->GetTouchableHandle() );
     G4VPhysicalVolume* oldTopVolume= fTrack->GetTouchableHandle()->GetVolume();
     G4VPhysicalVolume* newTopVolume=
     fNavigator->ResetHierarchyAndLocate( fTrack->GetPosition(), 
        fTrack->GetMomentumDirection(),
	*((G4TouchableHistory*)fTrack->GetTouchableHandle()()) );
//     if(newTopVolume != oldTopVolume ){
     if(newTopVolume != oldTopVolume || oldTopVolume->GetRegularStructureId() == 1 ) { 
        fTouchableHandle = fNavigator->CreateTouchableHistory();
        fTrack->SetTouchableHandle( fTouchableHandle );
        fTrack->SetNextTouchableHandle( fTouchableHandle );
     }
  }
// Set OriginTouchableHandle for primary track
   if(fTrack->GetParentID()==0){
     fTrack->SetOriginTouchableHandle(fTrack->GetTouchableHandle());
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
         G4cerr << "ERROR - G4SteppingManager::SetInitialStep()" << G4endl
                << "        Primary particle starting at - "
                << fTrack->GetPosition()
                << " - is outside of the world volume." << G4endl;
         G4Exception("G4SteppingManager::SetInitialStep()", "Tracking0010",
                     FatalException, "Primary vertex outside of the world!");
       }

       fTrack->SetTrackStatus( fStopAndKill );
       G4cout << "WARNING - G4SteppingManager::SetInitialStep()" << G4endl
              << "          Initial track position is outside world! - "
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

