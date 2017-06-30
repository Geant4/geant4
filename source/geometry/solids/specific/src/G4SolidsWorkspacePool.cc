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

#include "G4Types.hh"

#include "G4SolidsWorkspacePool.hh"
#include "G4SolidsWorkspace.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex singletonMutexSWP = G4MUTEX_INITIALIZER;
}

G4ThreadLocal G4SolidsWorkspace* G4SolidsWorkspacePool::fMyWorkspace=0;
G4SolidsWorkspacePool* G4SolidsWorkspacePool::thePool=0;

G4SolidsWorkspacePool* G4SolidsWorkspacePool::GetInstance()
{
  G4AutoLock l(&singletonMutexSWP);
  if( !thePool )  { thePool= new G4SolidsWorkspacePool(); }
  return thePool;
}

// For use with MT and current G4WorkerThread -- which uses static methods
//
G4SolidsWorkspace* G4SolidsWorkspacePool::CreateWorkspace()
{
  G4SolidsWorkspace* geometryWrk=0;
  if( !fMyWorkspace )
  {
    geometryWrk= new G4SolidsWorkspace();

    if( !geometryWrk )
    {
      G4Exception("GeometryWorspacePool::CreateWorkspace", "Geom-003",
                  FatalException, "Failed to create workspace.");
    }
    else
    {
      // geometryWrk->UseWorkspace();  // Do not assign it already.
      fMyWorkspace= geometryWrk;
    }
  }
  else
  {
    G4Exception("GeometryWorspacePool::CreateWorkspace", "Geom-003",
                FatalException,
                "Cannot create workspace twice for the same thread.");
    geometryWrk= fMyWorkspace; 
  }
  
  return geometryWrk;
}

// Create it (as above) and use it
//
void G4SolidsWorkspacePool::CreateAndUseWorkspace()
{
  (this->CreateWorkspace())->UseWorkspace();
}

// Reuse an existing workspace - or create a new one if needed.
//
G4SolidsWorkspace* G4SolidsWorkspacePool::FindOrCreateWorkspace()
{
  G4SolidsWorkspace* geometryWrk= fMyWorkspace;
  if( !geometryWrk )
  { 
    geometryWrk= this->CreateWorkspace(); 
  } 
  geometryWrk->UseWorkspace();
  
  fMyWorkspace= geometryWrk; // assign it for use by this thread.
  return geometryWrk;
}


// Keep the unused Workspace - for recycling
//
void G4SolidsWorkspacePool::Recycle( G4SolidsWorkspace *geometryWrk )
{
  geometryWrk->ReleaseWorkspace(); 
  delete geometryWrk;
}


void G4SolidsWorkspacePool::CleanUpAndDestroyAllWorkspaces()
{
  if (fMyWorkspace)
  {
     fMyWorkspace->DestroyWorkspace();
     delete fMyWorkspace;
     fMyWorkspace=0;
  }
}


G4SolidsWorkspacePool::G4SolidsWorkspacePool()
{
//  fWarehouse=0;
}


G4SolidsWorkspacePool::~G4SolidsWorkspacePool()
{
}
