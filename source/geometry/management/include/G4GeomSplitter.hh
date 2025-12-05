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
// G4GeomSplitter
//
// Class description:
//
// Utility template class for splitting of RW data for thread-safety from
// classes: G4LogicalVolume, G4Region, G4VPhysicalVolume, G4PolyconeSide
// G4PolyhedraSide, G4PVReplica. 

// Author: Xin Dong (Northeastern Univ.), 01.25.2009 - Initial version
// ------------------------------------------------------------------------
#ifndef G4GEOMSPLITTER_HH
#define G4GEOMSPLITTER_HH

#include "globals.hh"
#include "geomwdefs.hh"
#include "G4AutoLock.hh"

/**
 * @brief G4GeomSplitter is an utility class for splitting of R/W data
 * for thread-safety from geometry classes.
 * T is the private data from the object to be split
 */

template <class T>
class G4GeomSplitter
{
  public:

    /**
     * Constructor.
     */
    G4GeomSplitter()
      :  sharedOffset(nullptr)
    {
      G4MUTEXINIT(mutex);
    }

    /**
     * Reallocates data for a given 'size'.
     */
    T* Reallocate(G4int size)
    {
       totalspace = size;
       return (T *) std::realloc(offset, totalspace * sizeof(T));
    }

    /**
     * Invoked by the master or work thread to create a new subinstance
     * whenever a new split class instance is created.
     */
    G4int CreateSubInstance()
    {
      G4AutoLock l(&mutex);
      ++totalobj;
      if (totalobj > totalspace)
      {
        offset = Reallocate(totalspace+512);
        if (offset == nullptr)
        {
           G4Exception("G4GeomSPlitter::CreateSubInstance()",
                       "OutOfMemory", FatalException, "Cannot malloc space!");
        }
        sharedOffset = offset;
      }
      return (totalobj - 1);
    }

    /**
     * Utility to copy data from master in memory.
     */
    void CopyMasterContents()
    {
      G4AutoLock l(&mutex);
      std::memcpy(offset, sharedOffset, totalspace * sizeof(T));
    }
  
    /**
     * Invoked by each worker thread to copy all the subinstance array
     * from the master thread.
     */
    void SlaveCopySubInstanceArray()
    {
      G4AutoLock l(&mutex);
      if (offset != nullptr)  { return; }
      offset = Reallocate(totalspace);
      if (offset == nullptr)
      {
        G4Exception("G4GeomSplitter::SlaveCopySubInstanceArray()",
                    "OutOfMemory", FatalException, "Cannot malloc space!");
      }
      l.unlock();
      CopyMasterContents();
    }

    /**
     * Invoked by each worker thread to create the subinstance array and
     * initialize each subinstance using a particular method defined by
     * the subclass.
     */
    void SlaveInitializeSubInstance()
    {
      G4AutoLock l(&mutex);
      if (offset != nullptr)  { return; }
      offset = Reallocate(totalspace);

      if (offset == nullptr)
      {
        G4Exception("G4GeomSplitter::SlaveInitializeSubInstance()",
                    "OutOfMemory", FatalException, "Cannot malloc space!");
      }

      for (G4int i=0 ; i<totalspace; ++i)
      {
        offset[i].initialize();
      }
    }

    /**
     * Invoked by each worker thread at start of a run (2nd or later)
     * to copy again all the subinstance array from the master thread.
     * To cope with user's changes in Geometry - e.g. change of material
     * in a volume.
     */
    void SlaveReCopySubInstanceArray()
    {
      if (offset == nullptr)
      {
        SlaveInitializeSubInstance();
        G4Exception("G4GeomSPlitter::SlaveReCopySubInstance()",
                    "MissingInitialisation", JustWarning,
                    "Must be called after Initialisation or first Copy.");
      }
      CopyMasterContents();
    }
  
    /**
     * Invoked by all threads to free the subinstance array.
     */
    void FreeSlave()
    {
      if (offset == nullptr)  { return; }
      std::free( offset );
      offset = nullptr;
    }

    // Extension - to allow sharing of workspaces
  
    /**
     * Returns a pointer to the split data.
     */
    T* GetOffset() { return offset; }
  
    /**
     * Uses recycled work area - which was created previously.
     */
    void UseWorkArea( T* newOffset )
    {
      if( (offset!=nullptr) && (offset!=newOffset) )
      {
         G4Exception("G4GeomSplitter::UseWorkspace()", 
                     "TwoWorkspaces", FatalException,
                     "Thread already has workspace - cannot use another.");
      }
      offset = newOffset;
    }

    /**
     * Detaches the current thread from this location.
     * The object which calls this method is responsible for it.
     */
    T* FreeWorkArea()
    {
      T* offsetRet = offset;
      offset = nullptr;
      return offsetRet;
    }

  public:

    G4GEOM_DLL static G4ThreadLocal T* offset;

  private:

    G4int totalobj{0};
    G4int totalspace{0};
    T* sharedOffset;
    G4Mutex mutex;
};

template <typename T> G4ThreadLocal T* G4GeomSplitter<T>::offset = nullptr;

#endif
