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
// ------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
// Create and hold a pointer to Workspace.
// This class holds a thread-private static instance
// of the template parameter workspace.
//
// The concrete implementation of workspace objects
// are responsible for instantiating a singleton instance
// of this pool.
//
// Recycling of this pool can enable reuse among different
// threads in task-based - or 'on-demand' - simulation.

// ------------------------------------------------------------

#ifndef G4TWORKSPACEPOOL_HH
#define G4TWORKSPACEPOOL_HH

#include "tls.hh"
#include "globals.hh"

template<class T>
class G4TWorkspacePool
{
  public:

    inline T* CreateWorkspace();
    // For use with simple MT mode - each thread gets a workspace
    // and uses it until end
    
    inline void CreateAndUseWorkspace();
    // Create it (as above) and use it
    
    inline T* FindOrCreateWorkspace();
    // For use with 'dynamic' model of threading - workspaces can be recycled
    // Reuse an existing workspace - or create a new one if needed.
    // This will never fail, except if system is out of resources
    
    inline T* GetWorkspace() { return fMyWorkspace; }
    // Give back the existing, active workspace for my thread / task
    
    inline void Recycle( T * myWrkSpace );
    // Keep the unused Workspace - for recycling
    
    inline void CleanUpAndDestroyAllWorkspaces();
    // To be called once at the end of the job

  public:

    G4TWorkspacePool() {}
    ~G4TWorkspacePool() {}
    
  private:

    static G4ThreadLocal T* fMyWorkspace;
    // The thread's workspace - if assigned
};

template<typename T> G4ThreadLocal T* G4TWorkspacePool<T>::fMyWorkspace=0;

template<class T>
T* G4TWorkspacePool<T>::CreateWorkspace()
{
    T* wrk = 0;
    if ( !fMyWorkspace )
    {
        wrk = new T;
        if ( !wrk )
        {
            G4Exception("G4TWorspacePool<someType>::CreateWorkspace",
                        "MemoryError", FatalException,
                        "Failed to create workspace.");
        }
        else
        {
            fMyWorkspace = wrk;
        }
    }
    else
    {
        G4Exception("ParticlesWorspacePool::CreateWorkspace",
                    "InvalidCondition", FatalException,
                    "Cannot create workspace twice for the same thread.");
        wrk = fMyWorkspace;
    }
    return wrk;
}

template<class T>
void G4TWorkspacePool<T>::CreateAndUseWorkspace()
{
    (this->CreateWorkspace())->UseWorkspace();
}

template<class T>
T* G4TWorkspacePool<T>::FindOrCreateWorkspace()
{
    T* wrk= fMyWorkspace;
    if( !wrk )
    {
        wrk= this->CreateWorkspace();
    }
    wrk->UseWorkspace();
    
    fMyWorkspace= wrk; // assign it for use by this thread.
    return wrk;
}

template<class T>
void G4TWorkspacePool<T>::Recycle( T * myWrkSpace )
{
   myWrkSpace->ReleaseWorkspace(); 
   delete myWrkSpace;
}

template<class T>
void G4TWorkspacePool<T>::CleanUpAndDestroyAllWorkspaces()
{
   if (fMyWorkspace)
   {
      fMyWorkspace->DestroyWorkspace();
      delete fMyWorkspace;
      fMyWorkspace=0;
   }
}

#endif // G4TWORKSPACEPOOL_HH
