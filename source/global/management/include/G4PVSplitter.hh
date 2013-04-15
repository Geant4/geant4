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
// GEANT 4 class header file 
//
// Class Description:
//
// Utility template class for splitting of G4PhysicsVector data
// for thread-safety.

//      ---------------- G4PVSplitter ----------------
//
// Author: X.Dong (NorthEastern University), January 2009
// ------------------------------------------------------------

#ifndef G4PVSPLITTER_HH
#define G4PVSPLITTER_HH

#include <stdlib.h>

#include "G4Types.hh"

#ifdef G4MULTITHREADED
extern pthread_mutex_t mutexPhysicsVector;
// pthread_mutex_t mutexPhysicsVector = PTHREAD_MUTEX_INITIALIZER;
#endif

template <class T>  // T is the private data from the PV to be split
class G4PVSplitter
{
  public:
 
    G4PVSplitter() {}

    G4int CreateSubInstance()
      // Invoked by the master or work thread to create a new subinstance
      // whenever a new split class instance is created.
    {
#ifdef G4MULTITHREADED
      pthread_mutex_lock(&mutexPhysicsVector);
#endif

      totalobj++;
      if (totalobj > totalspace) NewSubInstances();

      if (phaseshadow == 0)
      {
        totalobjshadow = totalobj;
        offsetshadow = offset;
      }

      G4int totalobjlocal = totalobj;

#ifdef G4MULTITHREADED
      pthread_mutex_unlock(&mutexPhysicsVector);
#endif

      return (totalobjlocal - 1);
    }

    void SlaveInitializeSubInstance()
      // Invoked by each worker thread to create the subinstance array and
      // initialize each subinstance using a particular method defined by
      // the subclass.
    {
#ifdef G4MULTITHREADED
      pthread_mutex_lock(&mutexPhysicsVector);
#endif

      phaseshadow = 1;
      totalobj = totalobjshadow;
      totalspace = totalobj;
      offset = (T *) malloc(totalspace * sizeof(T));
      if (offset == 0)
      {
        G4Exception("G4PVSPlitter::SlaveInitializeSubInstance()",
                    "OutOfMemory", FatalException, "Cannot malloc space!");
      }

      // memcpy(offset, offsetshadow, totalspace * sizeof(T));

      for (G4int i = 0 ; i < totalspace ; i++)
      {
        offset[i].initialize();
      }

#ifdef G4MULTITHREADED
      pthread_mutex_unlock(&mutexPhysicsVector);
#endif
    }

    void NewSubInstances()
      // Invoked by each worker thread to grow the subinstance array and
      // initialize each new subinstance using a particular method defined
      // by the subclass.
    {
      if (totalspace >= totalobj)  { return; }
      G4int originaltotalspace = totalspace;
      totalspace = totalobj + 512;
      offset = (T *) realloc(offset, totalspace * sizeof(T));
      if (offset == 0)
      {
        G4Exception("G4PVSPlitter::NewSubInstances()", "OutOfMemory",
                    FatalException, "Cannot malloc space!");
          return;
      }

      for (G4int i = originaltotalspace ; i < totalspace ; i++)
      {
        offset[i].initialize();
      }
    }

    void FreeSlave()
      // Invoked by all threads to free the subinstance array.
    {
      if (!offset)  { return; }
      delete offset;
    }

  public:

    G4GLOB_DLL static T* offsetshadow;
    G4GLOB_DLL static G4ThreadLocal T* offset;

  private:

    G4GLOB_DLL static G4int phaseshadow; // 0: master init, 1: worker init
    G4GLOB_DLL static G4int totalobjshadow;
    G4GLOB_DLL static G4ThreadLocal G4int totalobj;
    G4GLOB_DLL static G4ThreadLocal G4int totalspace;
};

#endif
