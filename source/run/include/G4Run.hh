// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Run.hh,v 1.1 1999-01-07 16:14:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4Run_h
#define G4Run_h 1

#include "globals.hh"
#include "G4Allocator.hh"

class G4Event;
class G4HCtable;
class G4DCtable;

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

  public:
    inline void SetRunID(G4int id)
    { runID = id; }
    inline G4int GetRunID() const
    { return runID; }
    inline G4int GetNumberOfEvent() const
    { return numberOfEvent; }
    inline virtual void RecordEvent(G4Event*) 
    { numberOfEvent++; }
    inline void SetHCtable(G4HCtable* HCtbl)
    { HCtable = HCtbl; }
    inline const G4HCtable* GetHCtable() const
    { return HCtable; }
    inline void SetDCtable(G4DCtable* DCtbl)
    { DCtable = DCtbl; }
    inline const G4DCtable* GetDCtable() const
    { return DCtable; }
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

