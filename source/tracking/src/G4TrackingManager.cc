// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrackingManager.cc,v 1.6 2001-02-08 07:39:53 tsasaki Exp $
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
#include "G4Trajectory.hh"
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
  //  GimmeSecondaries()->clearAndDestroy();    
  //  G4std::vector<G4Track*>::iterator itr;
  G4int itr;
  //  for(itr=GimmeSecondaries()->begin();itr=GimmeSecondaries()->end();itr++){ 
  for(itr=0;itr<GimmeSecondaries()->size();itr++){ 
     delete GimmeSecondaries()->at(itr);
  }
  GimmeSecondaries()->clear();  
  
  // Pre tracking user intervention process.
  fpTrajectory = NULL;
  if( fpUserTrackingAction != NULL ) {
     fpUserTrackingAction->PreUserTrackingAction(fpTrack);
  }
#ifdef G4_STORE_TRAJECTORY
  // Construct a trajectory if it is requested
  if(StoreTrajectory&&(!fpTrajectory)) { 
     fpTrajectory = new G4Trajectory(fpTrack); // default trajectory concrete class object
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
     fpUserTrackingAction->PostUserTrackingAction(fpTrack);
  }

  // Destruct the trajectory if it was created
#ifdef G4VERBOSE
  if(StoreTrajectory&&verboseLevel>10) fpTrajectory->ShowTrajectory();
#endif
  if( (!StoreTrajectory)&&fpTrajectory ) {
      delete fpTrajectory;
  }
}

void G4TrackingManager::SetTrajectory(G4VTrajectory* aTrajectory)
{
#ifndef G4_STORE_TRAJECTORY
  G4Exception("G4TrackingManager::SetTrajectory is invoked without G4_STORE_TRAJECTORY compilor option");
#endif
  fpTrajectory = aTrajectory;
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
       G4cout << G4endl;
       G4cout << "*******************************************************"
            << "**************************************************"
            << G4endl;
       G4cout << "* G4Track Information: " 
            << "  Particle = " << fpTrack->GetDefinition()->GetParticleName() 
            << "," 
	    << "   Track ID = " << fpTrack->GetTrackID() 
            << "," 
	    << "   Parent ID = " << fpTrack->GetParentID() 
            << G4endl;
       G4cout << "*******************************************************"
            << "**************************************************"
            << G4endl;
       G4cout << G4endl;
      
     }
#endif
  }
}









