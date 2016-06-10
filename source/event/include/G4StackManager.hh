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
// $Id: G4StackManager.hh 66892 2013-01-17 10:57:59Z gunter $
//
//  Last Modification : 04/Oct/11 P. Mato - making use of G4TrackStack with value semantics
///


#ifndef G4StackManager_h
#define G4StackManager_h 1

#include "G4UserStackingAction.hh"
#include "G4StackedTrack.hh"
#include "G4TrackStack.hh"
#include "G4SmartTrackStack.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "globals.hh"
#include "evmandefs.hh"

class G4StackingMessenger;
class G4VTrajectory;

// class description:
//
// This is the manager class of handling stacks of G4Track objects.
// This class must be a singleton and be constructed by G4EventManager.
// Almost all methods must be invoked exclusively by G4EventManager.
// Especially, some Clear() methods MUST NOT be invoked by the user.
// Event abortion is handled by G4EventManager.
//
// This G4StackingManager has three stacks, the urgent stack, the
// waiting stack, and the postpone to next event stack. The meanings
// of each stack is descrived in the Geant4 user's manual.
//

class G4StackManager 
{
  public:
      G4StackManager();
      ~G4StackManager();

  private:
      const G4StackManager& operator=(const G4StackManager &right);
      G4int operator==(const G4StackManager &right) const;
      G4int operator!=(const G4StackManager &right) const;

  public:
      G4int PushOneTrack(G4Track *newTrack, G4VTrajectory *newTrajectory = 0);
      G4Track * PopNextTrack(G4VTrajectory**newTrajectory);
      G4int PrepareNewEvent();

  public: // with description
      void ReClassify();
      //  Send all tracks stored in the Urgent stack one by one to 
      // the user's concrete ClassifyNewTrack() method. This method
      // can be invoked from the user's G4UserStackingAction concrete
      // class, especially fron its NewStage() method. Be aware that
      // when the urgent stack becomes empty, all tracks in the waiting
      // stack are send to the urgent stack and then the user's NewStage()
      // method is invoked.

      void SetNumberOfAdditionalWaitingStacks(G4int iAdd);
      //  Set the number of additional (optional) waiting stacks.
      // This method must be invoked at PreInit, Init or Idle states.
      // Once the user set the number of additional waiting stacks,
      // he/she can use the corresponding ENUM in G4ClassificationOfNewTrack.
      // The user should invoke G4RunManager::SetNumberOfAdditionalWaitingStacks
      // method, which invokes this method.

      void TransferStackedTracks(G4ClassificationOfNewTrack origin, G4ClassificationOfNewTrack destination);
      //  Transfter all stacked tracks from the origin stack to the destination stack.
      // The destination stack needs not be empty.
      // If the destination is fKill, tracks are deleted.
      // If the origin is fKill, nothing happen.

      void TransferOneStackedTrack(G4ClassificationOfNewTrack origin, G4ClassificationOfNewTrack destination);
      //  Transfter one stacked track from the origin stack to the destination stack.
      // The transfered track is the one which came last to the origin stack.
      // The destination stack needs not be empty.
      // If the destination is fKill, the track is deleted.
      // If the origin is fKill, nothing happen.

  private:
      G4UserStackingAction * userStackingAction;
      G4int verboseLevel;
#ifdef G4_USESMARTSTACK
      G4SmartTrackStack * urgentStack;
#else
      G4TrackStack * urgentStack;
#endif
      G4TrackStack * waitingStack;
      G4TrackStack * postponeStack;
      G4StackingMessenger* theMessenger;
      std::vector<G4TrackStack*> additionalWaitingStacks;
      G4int numberOfAdditionalWaitingStacks;

  public:
      void clear();
      void ClearUrgentStack();
      void ClearWaitingStack(int i=0);
      void ClearPostponeStack();
      G4int GetNTotalTrack() const;
      G4int GetNUrgentTrack() const;
      G4int GetNWaitingTrack(int i=0) const;
      G4int GetNPostponedTrack() const;
      void SetVerboseLevel( G4int const value );
      void SetUserStackingAction(G4UserStackingAction* value);
  
  private:
     G4ClassificationOfNewTrack DefaultClassification(G4Track *aTrack);
};

#endif

