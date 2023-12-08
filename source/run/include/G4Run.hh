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
// G4Run
//
// Class description:
//
// This class represents a run. An object of this class is constructed
// and deleted by G4RunManager. Basically the user should use only the
// accessors (get methods). All properties are set by G4RunManager.

// Author: M.Asai, 1996
// --------------------------------------------------------------------
#ifndef G4Run_hh
#define G4Run_hh 1

#include "G4Profiler.hh"
#include "globals.hh"

#include <vector>

class G4Event;
class G4HCtable;
class G4DCtable;

class G4Run
{
  public:
    using ProfilerConfig = G4ProfilerConfig<G4ProfileType::Run>;

  public:
    G4Run();
    virtual ~G4Run();
    G4Run(const G4Run&) = delete;
    G4Run& operator=(const G4Run&) = delete;

    // Method to be overwritten by the user for recording events in this run.
    // In such a case, it is the user's responsibility to increment
    // numberOfEvent. Also, user's run class object must be instantiated in
    // user's runAction.
    virtual void RecordEvent(const G4Event*);

    // Method to be overwritten by the user for merging local G4Run object
    // to the global G4Run object.
    virtual void Merge(const G4Run*);

    // Method to be overwritten by the user for merging sub-event results
    // (e.g. hits) into the master G4Event.
    // [N.B.] Trajectories are merged in the base-class method. So, the user
    // should invoke the base-class method at the end of his/her overwriting
    // method.
    virtual void MergeSubEvent(G4Event* masterEv, const G4Event* subEv);

    // Store a G4Event object until this run object is deleted.
    // Given the potential large memory size of G4Event and its data-member
    // objects stored in G4Event, the user must be careful and responsible
    // for not storing too many G4Event objects. This method is invoked by
    // G4RunManager if the user invokes G4EventManager::KeepTheCurrentEvent()
    // or "/event/keepCurrentEvent" UI command while the particular event is
    // in being processed (typically in EndOfEventAction).
    void StoreEvent(G4Event* evt);

    // Returns the run ID. Run ID is set by G4RunManager.
    inline G4int GetRunID() const { return runID; }

    // Returns number of events processed in this run. The number is
    // incremented at the end of each event processing.
    inline G4int GetNumberOfEvent() const { return numberOfEvent; }

    inline G4int GetNumberOfEventToBeProcessed() const { return numberOfEventToBeProcessed; }

    // Returns hits collection.
    inline const G4HCtable* GetHCtable() const { return HCtable; }

    // Returns digi collection.
    inline const G4DCtable* GetDCtable() const { return DCtable; }

    // Returns random number status at the beginning of this run.
    inline const G4String& GetRandomNumberStatus() const { return randomNumberStatus; }

    // Returns the event vector.
    inline std::vector<const G4Event*>* GetEventVector() const { return eventVector; }

    inline void SetRunID(G4int id) { runID = id; }
    inline void SetNumberOfEventToBeProcessed(G4int n_ev) { numberOfEventToBeProcessed = n_ev; }
    inline void SetHCtable(G4HCtable* HCtbl) { HCtable = HCtbl; }
    inline void SetDCtable(G4DCtable* DCtbl) { DCtable = DCtbl; }
    inline void SetRandomNumberStatus(G4String& st) { randomNumberStatus = st; }

  protected:
    G4int runID = 0;
    G4int numberOfEvent = 0;
    G4int numberOfEventToBeProcessed = 0;
    G4HCtable* HCtable = nullptr;
    G4DCtable* DCtable = nullptr;
    G4String randomNumberStatus = "";
    std::vector<const G4Event*>* eventVector = nullptr;
};

#endif
