// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StackManager.hh,v 1.1 1999-01-07 16:06:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Last Modification : 09/Dec/96 M.Asai
//


#ifndef G4StackManager_h
#define G4StackManager_h 1

#include "G4UserStackingAction.hh"
#include "G4StackedTrack.hh"
#include "G4TrackStack.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "globals.hh"
class G4StackingMessenger;

class G4StackManager 
{
  public:
      G4StackManager();
      ~G4StackManager();

  private:
      const G4StackManager & operator=
                          (const G4StackManager &right);
      int operator==(const G4StackManager &right) const;
      int operator!=(const G4StackManager &right) const;

  public:
      G4int PushOneTrack(G4Track *newTrack);
      G4Track * PopNextTrack();
      G4int PrepareNewEvent();
      void ReClassify();

  private:
      G4UserStackingAction * userStackingAction;
      int verboseLevel;
      G4TrackStack * urgentStack;
      G4TrackStack * waitingStack;
      G4TrackStack * postponeStack;
      G4StackingMessenger* theMessenger;

  public:
      inline void clear()
      { 
        ClearUrgentStack();
        ClearWaitingStack();
      };
      inline void ClearUrgentStack()
      { urgentStack->clear(); };
      inline void ClearWaitingStack()
      { waitingStack->clear(); };
      inline void ClearPostponeStack()
      { postponeStack->clear(); };
      inline G4int GetNTotalTrack() const
      { return urgentStack->GetNTrack()
             + waitingStack->GetNTrack()
             + postponeStack->GetNTrack(); };
      inline G4int GetNUrgentTrack() const
      { return urgentStack->GetNTrack(); };
      inline G4int GetNWaitingTrack() const
      { return waitingStack->GetNTrack(); };
      inline G4int GetNPostponedTrack() const
      { return postponeStack->GetNTrack(); };
      inline void SetVerboseLevel( int const value )
      { verboseLevel = value; };
      inline void SetUserStackingAction(G4UserStackingAction* value)
      { 
        if (userStackingAction) delete userStackingAction;
	userStackingAction = value;
	userStackingAction->SetStackManager(this);
      };

  private:
      inline G4ClassificationOfNewTrack 
      DefaultClassification(G4Track *aTrack)
      { 
        G4ClassificationOfNewTrack classification = fUrgent;
        if( aTrack->GetTrackStatus() == fPostponeToNextEvent )
        { classification = fPostpone; }
        return classification;
      };
};

#endif

