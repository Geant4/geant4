// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Run.hh,v 1.3 1999-12-15 14:53:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4Run_h
#define G4Run_h 1

#include "globals.hh"
#include "G4Allocator.hh"

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
    inline void *operator new(size_t);
    inline void operator delete(void* aRun);

  protected:
    G4int runID;
    G4int numberOfEvent;
    G4HCtable* HCtable;
    G4DCtable* DCtable;

  public: // with description
    inline G4int GetRunID() const
    { return runID; }
    //  Returns the run ID. Run ID is set by G4RunManager.
    inline G4int GetNumberOfEvent() const
    { return numberOfEvent; }
    //  Returns number of events processed in this run. The number is
    // incremented at the end of each event processing.
    inline const G4HCtable* GetHCtable() const
    { return HCtable; }
    //  List of names of hits collection
    inline const G4DCtable* GetDCtable() const
    { return DCtable; }
    //  List of names of digi collection
  public:
    inline void SetRunID(G4int id)
    { runID = id; }
    inline virtual void RecordEvent(G4Event*) 
    { numberOfEvent++; }
    inline void SetHCtable(G4HCtable* HCtbl)
    { HCtable = HCtbl; }
    inline void SetDCtable(G4DCtable* DCtbl)
    { DCtable = DCtbl; }
};

extern G4Allocator<G4Run> aRunAllocator;

inline void* G4Run::operator new(size_t)
{
  void* aRun;
  aRun = (void*)aRunAllocator.MallocSingle();
  return aRun;
}

inline void G4Run::operator delete(void* aRun)
{
  aRunAllocator.FreeSingle((G4Run*)aRun);
}

#endif

