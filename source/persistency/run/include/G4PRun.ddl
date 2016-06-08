// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PRun.ddl,v 1.6 1999/12/15 14:51:27 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef G4PRun_h
#define G4PRun_h 1

#include "globals.hh"
#include "G4Allocator.hh"

#include "HepODBMS/odbms/HepODBMS.h"

#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

class G4Run;

class G4PRun
 : public HepPersObj
{
  public:
    G4PRun();
    G4PRun(const G4Run* aRun);
    virtual ~G4PRun();

    G4Run* MakeTransientObject();

  protected:
    G4Pint runID;
    G4Pint numberOfEvent;
//    G4HCtable* HCtable;
//    G4DCtable* DCtable;

  public:
    inline void SetRunID(G4int id)
    { runID = id; }
    inline G4int GetRunID() const
    { return runID; }
    inline G4int GetNumberOfEvent() const
    { return numberOfEvent; }
//    inline virtual void RecordEvent(G4Event*) 
//    { numberOfEvent++; }
//    inline void SetHCtable(G4HCtable* HCtbl)
//    { HCtable = HCtbl; }
//    inline const G4HCtable* GetHCtable() const
//    { return HCtable; }
//    inline void SetDCtable(G4DCtable* DCtbl)
//    { DCtable = DCtbl; }
//    inline const G4DCtable* GetDCtable() const
//    { return DCtable; }
};


#endif

