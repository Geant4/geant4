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
// $Id:$
//
// 
// class G4LogicalVolume
//
// Class description:
//
// Utility template class for splitting of RW data for thread-safety from
// classes: G4LogicalVolume, G4Region, G4VPhysicalVolume, G4PolyconeSide
// G4PolyhedraSide, G4PVReplica. 

// Author:
// 01.25.09 X.Dong: Initial version from automatic MT conversion.
// ------------------------------------------------------------------------
#ifndef G4GEOMSPLITTER_HH
#define G4GEOMSPLITTER_HH

#include <stdlib.h>

#include "globals.hh"
#include "geomwdefs.hh"
#include "G4AutoLock.hh"

template <class T>  // T is the private data from the object to be split
class G4GeomSplitter
{
  public:

    G4GeomSplitter()
      : totalobj(0), totalspace(0), sharedOffset(0)
    {
      G4MUTEXINIT(mutex);
    }

    G4int CreateSubInstance()
      // Invoked by the master or work thread to create a new subinstance
      // whenever a new split class instance is created.
    {
      G4AutoLock l(&mutex);
      totalobj++;
      if (totalobj > totalspace)
      {
        totalspace=totalspace + 512;
        offset = (T *) realloc(offset, totalspace * sizeof(T));
        if (offset == 0)
        {
           G4Exception("G4GeomSPlitter::CreateSubInstance()",
                       "OutOfMemory", FatalException, "Cannot malloc space!");
        }
        sharedOffset = offset;
      }
      return (totalobj - 1);
    }

    void CopyMasterContents()
    {
      G4AutoLock l(&mutex);
      memcpy(offset, sharedOffset, totalspace * sizeof(T));
    }
  
    void SlaveCopySubInstanceArray()
      // Invoked by each worker thread to copy all the subinstance array
      // from the master thread.
    {
      G4AutoLock l(&mutex);
      if (offset)  { return; }
      offset = (T *) realloc(offset, totalspace * sizeof(T));
      if (offset == 0)
      {
        G4Exception("G4GeomSplitter::SlaveCopySubInstanceArray()",
                    "OutOfMemory", FatalException, "Cannot malloc space!");
      }
      l.unlock();
      CopyMasterContents();
    }

    void SlaveInitializeSubInstance()
      // Invoked by each worker thread to create the subinstance array and
      // initialize each subinstance using a particular method defined by
      // the subclass.
    {
      G4AutoLock l(&mutex);
      if (offset)  { return; }
      offset = (T *) realloc(offset, totalspace * sizeof(T));

      if (offset == 0)
      {
        G4Exception("G4GeomSplitter::SlaveInitializeSubInstance()",
                    "OutOfMemory", FatalException, "Cannot malloc space!");
      }

      for (G4int i = 0 ; i < totalspace ; i++)
      {
        offset[i].initialize();
      }
    }

    void SlaveReCopySubInstanceArray()
      // Invoked by each worker thread at start of a run (2nd or later)
      // to copy again all the subinstance array from the master thread.
      // To cope with user's changes in Geometry - e.g. change of material
      // in a volume
    {
      if (!offset)
      {
        SlaveInitializeSubInstance();
        G4Exception("G4GeomSPlitter::SlaveReCopySubInstance()",
                    "MissingInitialisation", JustWarning,
                    "Must be called after Initialisation or first Copy.");
      }
      CopyMasterContents();
    }
  
    void FreeSlave()
      // Invoked by all threads to free the subinstance array.
    {
      if (!offset)  { return; }
      free( offset );
      offset = 0;
    }

    // Extension - to allow sharing of workspaces
  
    T* GetOffset() { return offset; }
  
    void UseWorkArea( T* newOffset )
      // Use recycled work area - which was created previously
    {
      if( offset && offset!=newOffset )
      {
         G4Exception("G4GeomSplitter::UseWorkspace()", 
                     "TwoWorkspaces", FatalException,
                     "Thread already has workspace - cannot use another.");
      }
      offset= newOffset;
    }

    T* FreeWorkArea()
      // Detach this thread from this Location.
      // The object which calls this method is responsible for it.
    {
      T* offsetRet= offset;
      offset= 0;
      return offsetRet;
    }

  public:

    G4GEOM_DLL static G4ThreadLocal T* offset;

  private:

    G4int totalobj;
    G4int totalspace;
    T* sharedOffset;
    G4Mutex mutex;
};

template <typename T> G4ThreadLocal T* G4GeomSplitter<T>::offset = 0;

#endif
