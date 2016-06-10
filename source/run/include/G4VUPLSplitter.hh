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
// $Id: $
//
// 
// ------------------------------------------------------------
//
//	GEANT 4 class header file 
//
// ---------------- G4UPLSplitter ----------------
//
// Utility template class for splitting RW data for thread-safety from
// classes: G4UserPhysicsList, G4VPhysicsConstructor and G4CModularPhsyicsList
//
// ------------------------------------------------------------
// History:
// 01.25.2009 Xin Dong: First implementation from automatic MT conversion.
// ------------------------------------------------------------
#ifndef G4VUPLSPLITTER_HH
#define G4VUPLSPLITTER_HH

#include <stdlib.h>

#include "globals.hh"
#include "rundefs.hh"

//
// This class implements the split-mechanism for shared objects.
// Let's see how it works.
// In the split-class we have an instance of this class and an G4int instanceID
// Every time, in the master thread a new instance of the split-class
// is created the constructor calls:
//   instanceID = g4vuplsplitter.CreateInstance();
// This creates in memory an "array", pointed by "sharedOffset" of capacity "totalspace"
// The array contains "totalobj" (<=totalspace) instances (i.e. the array has
// un-initialized spaces)
// Note that also the TLS variables "offset" and "slavetotalspace" have also the same stuff
// When a worker thread is started we can call g4vuplsplitter.NewSubInstances()
// This will simply allocate enough space in the TLS space "offset" and call
// T::initialize() onto the new created methods.
// Alternatively one can call, when the worker thread start, g4vuplsplitter.SlaveCopySubInstanceArray()
// That will copy the content of master thread "array" into the TLS one

// To see this stuff in action see:
// G4VUserPhysicsList class and G4WorkerThread classes.

template <class T>  // T is the private data from the object to be split
class G4VUPLSplitter
{
  public:

    G4VUPLSplitter() : totalobj(0),totalspace(0),sharedOffset(0) {}

    G4int CreateSubInstance()
      // Invoked by the master thread to create a new subinstance
      // whenever a new split class instance is created.
      // This is called by constructor of shared classes, thus only master thread
      // calls this
    {
        //One more instance
        totalobj++;
        //If the number of objects is larger than the available spaces,
        //a re-allocation is needed
        if (totalobj > slavetotalspace)  { NewSubInstances(); }
        //Since this is called by Master thread, we can remember this
        totalspace = slavetotalspace;
        sharedOffset = offset;
        return (totalobj - 1);
    }

    void NewSubInstances()
      // Invoked by each worker thread to grow the subinstance array and
      // initialize each new subinstance using a particular method defined
      // by the subclass.
    {
        if (slavetotalspace  >= totalobj)  { return; }
        //Remember current large size
        G4int originaltotalspace = slavetotalspace;
        //Increase its size by some value (purely arbitrary)
        slavetotalspace = totalobj + 512;
        //Now re-allocate new space
        offset = (T *) realloc(offset, slavetotalspace * sizeof(T));
        if (offset == 0)
        {
            G4Exception("G4VUPLSplitter::NewSubInstances()",
                        "OutOfMemory", FatalException, "Cannot malloc space!");
            return;
        }
        //The newly created objects need to be initialized
        for (G4int i = originaltotalspace; i < slavetotalspace; i++)
        {
            offset[i].initialize();
        }
    }

    void FreeSlave()
      // Invoked by all threads to free the subinstance array.
    {
      if (!offset)  { return; }
      free( offset);
      offset = 0;
    }

    void SlaveCopySubInstanceArray()
    //Invoked by each worker thread to copy all subinstances array from
    //the master thread
    {
        if ( offset ) return;
        //Since this is called by worker threds, totalspace is some valid number > 0
        //Remember totalspace is the number of availabel slots from master.
        //We are sure that it has valid data
        offset = (T *)realloc(offset,totalspace * sizeof(T));
        if (offset == 0)
        {
            G4Exception("G4VUPLSplitter::SlaveCopySubInstanceArray()",
                        "OutOfMemory", FatalException, "Cannot malloc space!");
            return;
        }
        //Now just copy from master thread (sharedOffset)
        memcpy(offset,sharedOffset,totalspace*sizeof(T));
    }
  public:

    G4RUN_DLL static G4ThreadLocal G4int slavetotalspace; //Per-thread available number of slots
    G4RUN_DLL static G4ThreadLocal T* offset; //Pointer to first instance of an array

  private:

    G4int totalobj; //Total number of instances from master thread
    G4int totalspace; // Available number of "slots"
    T* sharedOffset;
};

#endif
