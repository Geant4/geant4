// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrackingManager.cc,v 1.1 1999-01-07 16:14:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4TrackingManager.cc
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#include "G4TrackingManager.hh"
#include "G4TrackingMessenger.hh"
#include "G4ios.hh"

G4TrackingMessenger * messenger;
//////////////////////////////////////
G4TrackingManager::G4TrackingManager()
//////////////////////////////////////
  : verboseLevel(0),StoreTrajectory(false),
    fpTrajectory(NULL), fpUserTrackingAction(NULL) 
{
  fpSteppingManager = new G4SteppingManager();
  messenger = new G4TrackingMessenger(this);
}

///////////////////////////////////////
G4TrackingManager::~G4TrackingManager()
///////////////////////////////////////
{
  delete messenger;
  delete fpSteppingManager;
  if (fpUserTrackingAction) delete fpUserTrackingAction;
}

////////////////////////////////////////////////////////////////
void G4TrackingManager::ProcessOneTrack(G4Track* apValueG4Track)
////////////////////////////////////////////////////////////////
{

  // Receiving a G4Track from the EventManager, this funciton has the
  // responsibility to trace the track till it stops.
  fpTrack = apValueG4Track;

  // Clear 2ndary particle vector
  GimmeSecondaries()->clearAndDestroy();    

  // Pre tracking user intervention process.
  if( fpUserTrackingAction != NULL ) {
     fpUserTrackingAction->PreUserTrackingAction();
  }
#ifdef G4_STORE_TRAJECTORY
  // Construct a trajectory if it is requested
  if(StoreTrajectory) { 
     fpTrajectory = new G4Trajectory(fpTrack); 
  }
  else { 
     fpTrajectory = NULL; 
  }
#endif

#ifdef G4VERBOSE
                         // !!!!! Verbose
                         if(verboseLevel>0) Verbose("ProcessOneTrack");
#endif

  // Give SteppingManger the pointer to the track which will be tracked 
  fpSteppingManager->SetInitialStep(fpTrack);

  // Give SteppingManger the maxmimum number of processes 
  fpSteppingManager->GetProcessNumber();

  // Give track the pointer to the Step
  fpTrack->SetStep(fpSteppingManager->GetStep());

  // Inform beginning of tracking to physics processes 
  fpTrack->GetDefinition()->GetProcessManager()->StartTracking();

  // Track the particle Step-by-Step while it is alive
  G4StepStatus stepStatus;

  while( (fpTrack->GetTrackStatus() == fAlive) ||
         (fpTrack->GetTrackStatus() == fStopButAlive) ){

    fpTrack->IncrementCurrentStepNumber();
    stepStatus = fpSteppingManager->Stepping();
#ifdef G4_STORE_TRAJECTORY
    if(StoreTrajectory) fpTrajectory->
                        AppendStep(fpSteppingManager->GetStep()); 
#endif
  }
  // Inform end of tracking to physics processes 
  fpTrack->GetDefinition()->GetProcessManager()->EndTracking();

  // Post tracking user intervention process.
  if( fpUserTrackingAction != NULL ) {
     fpUserTrackingAction->PostUserTrackingAction();
  }

  // Destruct the trajectory if it was created
#ifdef G4VERBOSE
  if(StoreTrajectory&&verboseLevel>10) fpTrajectory->ShowTrajectory();
#endif
  if( (!StoreTrajectory)&&fpTrajectory ) {
      delete fpTrajectory;
  }
}

//////////////////////////////////////
void G4TrackingManager::EventAborted()
//////////////////////////////////////
{
}

//************************************************************************
//
//  Private Member Functions
//
//************************************************************************


////////////////////////////////////////////////
void G4TrackingManager::Verbose(G4String select)
////////////////////////////////////////////////
{

  // !!!!! Verbose
  if( select == "ProcessOneTrack" ){
#ifdef G4VERBOSE
     if(verboseLevel >= 1) { 
       G4cout << endl;
       G4cout << "*******************************************************"
            << "**************************************************"
            << endl;
       G4cout << "* G4Track Information: " 
            << "  Particle = " << fpTrack->GetDefinition()->GetParticleName() 
            << "," 
	    << "   Track ID = " << fpTrack->GetTrackID() 
            << "," 
	    << "   Parent ID = " << fpTrack->GetParentID() 
            << endl;
       G4cout << "*******************************************************"
            << "**************************************************"
            << endl;
       G4cout << endl;
      
     }
#endif
  }
}









