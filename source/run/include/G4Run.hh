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
// $Id: G4Run.hh 68798 2013-04-05 13:28:19Z gcosmo $
//

#ifndef G4Run_h
#define G4Run_h 1

#include "globals.hh"
#include <vector>
class G4Event;
class G4HCtable;
class G4DCtable;

// class description:
//
//  This class represents a run. An object of this class is constructed
// and deleted by G4RunManager. Basically the user should use only the
// get methods. All properties are set by G4RunManager.
//  

class G4Run
{
  public:
    G4Run();
    virtual ~G4Run();

  private:
    // These copy constructor and = operator must not be used.
    G4Run(const G4Run &) {;}
    G4Run& operator=(const G4Run &) { return *this; }

  protected:
    G4int runID;
    G4int numberOfEvent;
    G4int numberOfEventToBeProcessed;
    G4HCtable* HCtable;
    G4DCtable* DCtable;
    G4String randomNumberStatus;
    std::vector<const G4Event*>* eventVector;

  public: // with description
    virtual void RecordEvent(const G4Event*);
    //  Method to be overwritten by the user for recording events in this run.
    //  In such a case, it is the user's responsibility to increment numberOfEvent.
    //  Also, user's run class object must be instantiated in user's runAction.
    virtual void Merge(const G4Run*);
    //  Method to be overwritten by the user for merging local G4Run object to 
    //  the global G4Run object.

  public: // with description
    inline G4int GetRunID() const
    { return runID; }
    //  Returns the run ID. Run ID is set by G4RunManager.
    inline G4int GetNumberOfEvent() const
    { return numberOfEvent; }
    //  Returns number of events processed in this run. The number is
    // incremented at the end of each event processing.
    inline G4int GetNumberOfEventToBeProcessed() const
    { return numberOfEventToBeProcessed; }
    inline const G4HCtable* GetHCtable() const
    { return HCtable; }
    //  List of names of hits collection
    inline const G4DCtable* GetDCtable() const
    { return DCtable; }
    //  List of names of digi collection
    inline const G4String& GetRandomNumberStatus() const
    { return randomNumberStatus; }
    // Return random number status at the beginning of this run
  public:
    inline void SetRunID(G4int id)
    { runID = id; }
    inline void SetNumberOfEventToBeProcessed(G4int n_ev)
    { numberOfEventToBeProcessed = n_ev; }
    inline void SetHCtable(G4HCtable* HCtbl)
    { HCtable = HCtbl; }
    inline void SetDCtable(G4DCtable* DCtbl)
    { DCtable = DCtbl; }
    inline void SetRandomNumberStatus(G4String& st)
    { randomNumberStatus = st; }

  public: // with description
    void StoreEvent(G4Event* evt);
    // Store a G4Event object until this run object is deleted.
    // Given the potential large memory size of G4Event and its datamember
    // objects stored in G4Event, the user must be careful and responsible for
    // not to store too many G4Event objects. This method is invoked by G4RunManager
    // if the user invokes G4EventManager::KeepTheCurrentEvent() or 
    // /event/keepCurrentEvent UI command while the particular event is in process
    // (typically in EndOfEventAction).
    inline const std::vector<const G4Event*>* GetEventVector() const
    { return eventVector; }
    // Return the event vector
};


#endif

