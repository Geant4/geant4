// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserStackingAction.hh,v 1.1 1999-01-07 16:06:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Last Modification : 09/Dec/96 M.Asai
//

#ifndef G4UserStackingAction_h
#define G4UserStackingAction_h 1

class G4StackManager;
class G4Track;
#include "G4ClassificationOfNewTrack.hh"

class G4UserStackingAction 
{
  public:
      G4UserStackingAction();
      virtual ~G4UserStackingAction();
  protected:
      G4StackManager * stackManager;
  public:
      inline void SetStackManager(G4StackManager * value)
      { stackManager = value; };

  public:
//---------------------------------------------------------------
// vitual methods to be implemented by user
//---------------------------------------------------------------
//
      virtual G4ClassificationOfNewTrack 
        ClassifyNewTrack(G4Track *const aTrack);
//
//    Reply G4ClassificationOfNewTrack determined by the
//  newly coming G4Track.
//
//    enum G4ClassificationOfNewTrack
//    {
//      fUrgent,    // put into the urgent stack
//      fWaiting,   // put into the waiting stack
//      fPostpone,  // postpone to the next event
//      fKill       // kill without stacking
//    };
//
//    The parent_ID of the track indicates the origin of it.
//                
//    G4int parent_ID = aTrack->get_parentID();
//   
//      parent_ID = 0 : primary particle
//                > 0 : secondary particle
//                < 0 : postponed from the previous event
//
//---------------------------------------------------------------
//
      virtual void NewStage();
//
//    This method is called by G4StackManager when the urgentStack
//  becomes empty and contents in the waitingStack are transtered
//  to the urgentStack.
//    Note that this method is not called at the begining of each
//  event, but "PrepareNewEvent" is called.
//
//    In case re-classification of the stacked tracks is needed,
//  use the following method to request to G4StackManager.
//
//    stackManager->ReClassify();
//
//  All of the stacked tracks in the waitingStack will be re-classified 
//  by "ClassifyNewTrack" method.
//    To abort current event, use the following method.
//
//    stackManager->clear();
//
//  Note that this way is valid and safe only for the case it is called
//  from this user class. The more global way of event abortion is
//
//    G4UImanager * UImanager = G4UImanager::GetUIpointer();
//    UImanager->ApplyCommand("/event/abort");
//
//---------------------------------------------------------------
//
      virtual void PrepareNewEvent();
//
//    This method is called by G4StackManager at the begining of
//  each event.
//    Be careful that the urgentStack and the waitingStack of 
//  G4StackManager are empty at this moment, because this method
//  is called before accepting primary particles. Also, note that
//  the postponeStack of G4StackManager may have some postponed
//  tracks.
//
//---------------------------------------------------------------

};

#endif

