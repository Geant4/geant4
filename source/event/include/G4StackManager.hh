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
// G4StackManager
//
// Class description:
//
// This is the manager class of handling stacks of G4Track objects.
// This class must be a singleton and be constructed by G4EventManager.
// Almost all methods must be invoked exclusively by G4EventManager.
// Especially, some Clear() methods MUST NOT be invoked by the user.
// Event abortion is handled by G4EventManager.
//
// G4StackManager has three stacks, the urgent stack, the
// waiting stack, and the postpone to next event stack. The meanings
// of each stack is descrived in the Geant4 User's Manual.

// Author: Makoto Asai, 1996
//
// History:
// - 01/Feb/1996, Makoto Asai - Created
// - 04/Oct/2011, Pere Mato - Use of G4TrackStack with value semantics
// - 28/Aug/2023, Makoto Asai - Adding sub-event parallelism
// --------------------------------------------------------------------
#ifndef G4StackManager_hh
#define G4StackManager_hh 1

#include <map>
#include <vector>

#include "G4UserStackingAction.hh"
#include "G4StackedTrack.hh"
#include "G4TrackStack.hh"
#include "G4SmartTrackStack.hh"
#include "G4SubEventTrackStack.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ExceptionSeverity.hh"
#include "globals.hh"

class G4StackingMessenger;
class G4VTrajectory;
class G4Event;
class G4ParticleDefinition;

class G4StackManager 
{
  public:

    G4StackManager();
   ~G4StackManager();

    const G4StackManager& operator=(const G4StackManager&) = delete;
    G4bool operator==(const G4StackManager&) const = delete;
    G4bool operator!=(const G4StackManager&) const = delete;

    G4int PushOneTrack(G4Track* newTrack,
                       G4VTrajectory* newTrajectory = nullptr);
    G4Track* PopNextTrack(G4VTrajectory** newTrajectory);
    G4int PrepareNewEvent(G4Event* currentEvent);

    void ReClassify();
      // Send all tracks stored in the Urgent stack one by one to 
      // the user's concrete ClassifyNewTrack() method. This method
      // can be invoked from the user's G4UserStackingAction concrete
      // class, especially fron its NewStage() method. Be aware that
      // when the urgent stack becomes empty, all tracks in the waiting
      // stack are send to the urgent stack and then the user's NewStage()
      // method is invoked.

    void SetNumberOfAdditionalWaitingStacks(G4int iAdd);
      // Set the number of additional (optional) waiting stacks.
      // This method must be invoked at PreInit, Init or Idle states.
      // Once the user set the number of additional waiting stacks,
      // he/she can use the corresponding ENUM in G4ClassificationOfNewTrack.
      // The user should invoke G4RunManager::SetNumberOfAdditionalWaitingStacks
      // method, which invokes this method.

    void TransferStackedTracks(G4ClassificationOfNewTrack origin,
                               G4ClassificationOfNewTrack destination);
      // Transfer all stacked tracks from the origin stack to the
      // destination stack. The destination stack needs not be empty.
      // If the destination is fKill, tracks are deleted.
      // If the origin is fKill, nothing happen.

    void TransferOneStackedTrack(G4ClassificationOfNewTrack origin,
                                 G4ClassificationOfNewTrack destination);
      // Transfter one stacked track from the origin stack to the destination
      // stack.
      // The transfered track is the one which came last to the origin stack.
      // The destination stack needs not be empty.
      // If the destination is fKill, the track is deleted.
      // If the origin is fKill, nothing happen.

    void RegisterSubEventType(G4int ty, G4int maxEnt);
      // Registering a sub-event type and the capacity of the tracks to be 
      // stored in a G4SubEvent object of the corresponding sub-event type

    void SetDefaultClassification(
             G4TrackStatus, G4ClassificationOfNewTrack,
             G4ExceptionSeverity es = G4ExceptionSeverity::IgnoreTheIssue);
    void SetDefaultClassification(
             const G4ParticleDefinition*, G4ClassificationOfNewTrack,
             G4ExceptionSeverity es = G4ExceptionSeverity::IgnoreTheIssue);
      // Define the default classification for a newly arriving track.
      // Default can be alternated by the UserStackingAction.
      // G4ExceptionSeverity can be set to warn the user if the classification
      // is inproperly changed.

    inline G4ClassificationOfNewTrack GetDefaultClassification()
    { return fDefaultClassification; }

  public:
    void ReleaseSubEvent(G4int ty);
    inline std::size_t GetNSubEventTypes()
    { return subEvtTypes.size(); }
    inline G4int GetSubEventType(std::size_t i)
    { return subEvtTypes[i]; }
    
    void clear();
    void ClearUrgentStack();
    void ClearWaitingStack(G4int i=0);
    void ClearPostponeStack();
    G4int GetNTotalTrack() const;
    G4int GetNUrgentTrack() const;
    G4int GetNWaitingTrack(G4int i=0) const;
    G4int GetNPostponedTrack() const;
    void SetVerboseLevel( G4int const value );
    void SetUserStackingAction(G4UserStackingAction* value);

  private:
    void DefineDefaultClassification(const G4Track* aTrack);
    void SortOut(G4StackedTrack&,G4ClassificationOfNewTrack);

  private:

    G4UserStackingAction* userStackingAction = nullptr;
    G4int verboseLevel = 0;
#ifdef G4_USESMARTSTACK
    G4SmartTrackStack* urgentStack = nullptr;
#else
    G4TrackStack* urgentStack = nullptr;
#endif
    G4TrackStack* waitingStack = nullptr;
    G4TrackStack* postponeStack = nullptr;
    G4StackingMessenger* theMessenger = nullptr;
    std::vector<G4TrackStack*> additionalWaitingStacks;
    G4int numberOfAdditionalWaitingStacks = 0;

    std::map<G4TrackStatus,
             std::pair<G4ClassificationOfNewTrack,G4ExceptionSeverity>>
             defClassTrackStatus;
    std::map<const G4ParticleDefinition*,
             std::pair<G4ClassificationOfNewTrack,G4ExceptionSeverity>>
             defClassPartDef;
    G4ClassificationOfNewTrack fDefaultClassification = fUrgent;
    G4ExceptionSeverity fExceptionSeverity = G4ExceptionSeverity::IgnoreTheIssue;

    std::map<G4int,G4SubEventTrackStack*> subEvtStackMap;
    std::vector<G4int> subEvtTypes;
};

#endif
