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
// $Id: G4PRun.ddl,v 1.9 2001/07/11 10:02:26 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

// Class Description:
//
//      This is a class which represent a persistent run in Geant4.
//    A persistent run object is constructed by
//    G4PersistentRunMan::Store() or read back from a database
//    by G4PersistentRunMan::Retrieve().
//      G4PersistentRunMan::Store() passes a pointer of a transient
//    G4Run object.  The G4PRun object construct itself by copying
//    the data member of G4Run.
//

#ifndef G4PRun_h
#define G4PRun_h 1

#include "G4Pglobals.hh"
#include "G4Allocator.hh"

#include "HepODBMS/odbms/HepODBMS.h"

#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

class G4Run;

class G4PRun
 : public HepPersObj
{
  public: // with description
    G4PRun();
    G4PRun(const G4Run* aRun);
    // Constructors.
    virtual ~G4PRun();
    // Destructor.
    G4Run* MakeTransientObject();
    // Construct a transient G4Run object from the data members of
    // G4PRun.

  protected:
    G4Pint runID;
    G4Pint numberOfEvent;
//    G4HCtable* HCtable;
//    G4DCtable* DCtable;

  public: // with description
    inline void SetRunID(G4int id)
    { runID = id; }
    inline G4int GetRunID() const
    { return runID; }
    // set and get Run ID.
    inline G4int GetNumberOfEvent() const
    { return numberOfEvent; }
    // get number of events in this run.

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

