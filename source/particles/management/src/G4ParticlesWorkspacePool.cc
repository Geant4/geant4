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

#include "G4ParticlesWorkspacePool.hh"
#include "G4ParticlesWorkspace.hh"

#include "G4AutoLock.hh"
namespace {
    G4Mutex singletonMutexSWP = G4MUTEX_INITIALIZER;
}

G4ThreadLocal G4ParticlesWorkspace* G4ParticlesWorkspacePool::fMyWorkspace=0;

G4ParticlesWorkspacePool* G4ParticlesWorkspacePool::thePool=0;

// static
G4ParticlesWorkspacePool* G4ParticlesWorkspacePool::GetInstance()
{
   G4AutoLock l(&singletonMutexSWP);
   if( !thePool ) thePool= new G4ParticlesWorkspacePool();
   return thePool;
}

// For use with MT and current G4WorkerThread -- which uses static methods
G4ParticlesWorkspace* G4ParticlesWorkspacePool::CreateWorkspace()
{
  G4ParticlesWorkspace* particlesWrk=0;
  if( !fMyWorkspace ){
    particlesWrk= new G4ParticlesWorkspace();

    if( !particlesWrk ) {
      G4Exception("ParticlesWorspacePool::CreateWorkspace", "Run0035",
                  FatalException, "Failed to create workspace.");
    }else{
       fMyWorkspace= particlesWrk;
    }
  }else{
    G4Exception("ParticlesWorspacePool::CreateWorkspace", "Run0035",
                FatalException,
                "Cannot create workspace twice for the same thread.");
    particlesWrk= fMyWorkspace; 
  }
  
  return particlesWrk;
}


void G4ParticlesWorkspacePool::CreateAndUseWorkspace()  // Create it (as above) and use it
{
  (this->CreateWorkspace())->UseWorkspace();
}

// Reuse an existing workspace - or create a new one if needed.
//
G4ParticlesWorkspace* G4ParticlesWorkspacePool::FindOrCreateWorkspace()
{
   G4ParticlesWorkspace* particlesWrk= fMyWorkspace;
   if( !particlesWrk ){ 
      particlesWrk= this->CreateWorkspace(); 
   } 
   particlesWrk->UseWorkspace();
  
   fMyWorkspace= particlesWrk; // assign it for use by this thread.
   return particlesWrk;
}

#if 0
void G4ParticlesWorkspacePool::ReleaseAndDestroyMyWorkspace()
{
  ReleaseAndDestroyWorkspace(fMyWorkspace);
  fMyWorkspace=0;
}

void G4ParticlesWorkspacePool::ReleaseAndDestroyWorkspace(G4ParticlesWorkspace *particlesWrk)
{
   particlesWrk->ReleaseWorkspace(); 
   delete particlesWrk; 
}
#endif 

void G4ParticlesWorkspacePool::Recycle( G4ParticlesWorkspace *particlesWrk )
{
   particlesWrk->ReleaseWorkspace(); 
//   if( fWarehouse ){ 
//   } else {
    delete particlesWrk;
//   }
}
      // Keep the unused Workspace - for recycling


void G4ParticlesWorkspacePool::CleanUpAndDestroyAllWorkspaces()
{
}


G4ParticlesWorkspacePool::G4ParticlesWorkspacePool()
{
//  fWarehouse=0;
}

G4ParticlesWorkspacePool::~G4ParticlesWorkspacePool()
{
}



