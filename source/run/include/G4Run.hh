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
// $Id: G4Run.hh,v 1.11 2005/09/19 16:53:53 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#ifndef G4Run_h
#define G4Run_h 1

#include "globals.hh"

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

  protected:
    G4int runID;
    G4int numberOfEvent;
    G4int numberOfEventToBeProcessed;
    G4HCtable* HCtable;
    G4DCtable* DCtable;

  public: // with description
    virtual void RecordEvent(const G4Event*);
    //  Method to be overwritten by the user for recording events in this run.
    //  In such a case, it is the user's responsibility to increment numberOfEvent.
    //  Also, user's run class object must be instantiated in user's runAction.

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
  public:
    inline void SetRunID(G4int id)
    { runID = id; }
    inline void SetNumberOfEventToBeProcessed(G4int n_ev)
    { numberOfEventToBeProcessed = n_ev; }
    inline void SetHCtable(G4HCtable* HCtbl)
    { HCtable = HCtbl; }
    inline void SetDCtable(G4DCtable* DCtbl)
    { DCtable = DCtbl; }
};


#endif

