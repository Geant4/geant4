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
// G4VUPLSplitter
//
// Class description:
//
// Utility template class for splitting RW data for thread-safety from classes:
// G4UserPhysicsList, G4VPhysicsConstructor and G4VModularPhysicsList.
// This class implements the split-mechanism for shared objects.
// In the split-class we have an instance of this class and an 'instanceID'.
// Every time in the master thread a new instance of the split-class is
// created, the constructor calls:
//   instanceID = g4vuplsplitter.CreateInstance();
// This creates in memory an "array", pointed by "sharedOffset" of capacity
// "totalspace". The array contains "totalobj" (<=totalspace) instances
// (i.e. the array has un-initialized spaces). Note that also the TLS variables
// "offset" and "workertotalspace" have also the same stuff. When a worker
// thread is started we can call g4vuplsplitter.NewSubInstances(). This will
// simply allocate enough space in the TLS space "offset" and call
// T::initialize() onto the new created methods. Alternatively one can call,
// when the worker thread start, g4vuplsplitter.workerCopySubInstanceArray(),
// that will copy the content of master thread "array" into the TLS one.
// To see this stuff in action see the G4VUserPhysicsList and G4WorkerThread
// classes.

// Author: Xin Dong, 25 January 2009 - First implementation from
//                                     automatic MT conversion.
// --------------------------------------------------------------------
#ifndef G4VUPLSplitter_hh
#define G4VUPLSplitter_hh 1

#include "G4AutoLock.hh"
#include "globals.hh"

#include "rundefs.hh"
#include <stdlib.h>

template<class T>  // T is the private data from the object to be split
class G4VUPLSplitter
{
  public:
    G4VUPLSplitter() { G4MUTEXINIT(mutex); }

    // Invoked by the master thread to create a new subinstance
    // whenever a new split class instance is created.
    // This is called by constructor of shared classes,
    // thus only master thread calls this
    G4int CreateSubInstance()
    {
      G4AutoLock l(&mutex);
      // One more instance
      ++totalobj;
      // If the number of objects is larger than the available spaces,
      // a re-allocation is needed
      if (totalobj > workertotalspace) {
        l.unlock();
        NewSubInstances();
        l.lock();
      }
      // Since this is called by Master thread, we can remember this
      totalspace = workertotalspace;
      sharedOffset = offset;
      return (totalobj - 1);
    }

    // Invoked by each worker thread to grow the subinstance array and
    // initialize each new subinstance using a particular method defined
    // by the subclass.
    void NewSubInstances()
    {
      G4AutoLock l(&mutex);
      if (workertotalspace >= totalobj) {
        return;
      }
      // Remember current large size
      G4int originaltotalspace = workertotalspace;
      // Increase its size by some value (purely arbitrary)
      workertotalspace = totalobj + 512;
      // Now re-allocate new space
      offset = (T*)realloc(offset, workertotalspace * sizeof(T));
      if (offset == nullptr) {
        G4Exception("G4VUPLSplitter::NewSubInstances()", "OutOfMemory", FatalException,
                    "Cannot malloc space!");
        return;
      }
      // The newly created objects need to be initialized
      for (G4int i = originaltotalspace; i < workertotalspace; ++i) {
        offset[i].initialize();
      }
    }

    // Invoked by all threads to free the subinstance array.
    void FreeWorker()
    {
      if (offset == nullptr) {
        return;
      }
      free(offset);
      offset = nullptr;
    }

    T* GetOffset() { return offset; }

    void UseWorkArea(T* newOffset)
    {
      // Use recycled work area - which was created previously
      if (offset != nullptr && offset != newOffset) {
        G4Exception("G4VUPLSplitter::UseWorkspace()", "TwoWorkspaces", FatalException,
                    "Thread already has workspace - cannot use another.");
      }
      offset = newOffset;
    }

    T* FreeWorkArea()
    {
      // Detach this thread from this Location
      // The object which calls this method is responsible for it.
      //
      T* offsetRet = offset;
      offset = nullptr;

      return offsetRet;
    }

    // Invoked by each worker thread to copy all subinstances array from
    // the master thread
    void WorkerCopySubInstanceArray()
    {
      if (offset != nullptr) return;

      // Since this is called by worker threds, totalspace is some valid
      // number > 0. Remember totalspace is the number of available slots
      // from master. We are sure that it has valid data
      G4AutoLock l(&mutex);
      offset = (T*)realloc(offset, totalspace * sizeof(T));
      if (offset == nullptr) {
        G4Exception("G4VUPLSplitter::WorkerCopySubInstanceArray()", "OutOfMemory", FatalException,
                    "Cannot malloc space!");
        return;
      }
      // Now just copy from master thread (sharedOffset)
      std::memcpy(offset, sharedOffset, totalspace * sizeof(T));
    }

  public:
    // Per-thread available number of slots
    G4RUN_DLL G4ThreadLocalStatic G4int workertotalspace;
    // Pointer to first instance of an array
    G4RUN_DLL G4ThreadLocalStatic T* offset;

  private:
    G4int totalobj = 0;  // Total number of instances from master thread
    G4int totalspace = 0;  // Available number of "slots"
    T* sharedOffset = nullptr;
    G4Mutex mutex;
};

template<typename T>
G4ThreadLocal G4int G4VUPLSplitter<T>::workertotalspace = 0;
template<typename T>
G4ThreadLocal T* G4VUPLSplitter<T>::offset = nullptr;

#endif
