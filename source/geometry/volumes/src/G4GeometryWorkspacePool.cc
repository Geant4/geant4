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
// $Id: G4GeometryWorkspacePool.cc 103041 2017-03-10 11:47:01Z gcosmo $
//
// 
// Class G4GeometryWorkspacePool - implementation
//
// ----------------------------------------------------------------------

#include "G4GeometryWorkspacePool.hh"
#include "G4GeometryWorkspace.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex singletonM = G4MUTEX_INITIALIZER;
}

G4ThreadLocal G4GeometryWorkspace* G4GeometryWorkspacePool::fMyWorkspace=0;

G4GeometryWorkspacePool* G4GeometryWorkspacePool::thePool=0;

// ----------------------------------------------------------------------
//
G4GeometryWorkspacePool* G4GeometryWorkspacePool::GetInstance()
{
  G4AutoLock l(&singletonM);
  if ( !thePool )  { thePool= new G4GeometryWorkspacePool(); }
  return thePool;
}

// ----------------------------------------------------------------------
// For use with MT and current G4WorkerThread -- which uses static methods
//
G4GeometryWorkspace* G4GeometryWorkspacePool::CreateWorkspace()
{
  G4GeometryWorkspace* geometryWrk=0;
  if( !fMyWorkspace )
  {
    geometryWrk= new G4GeometryWorkspace();

    if( !geometryWrk )
    {
      G4Exception("GeometryWorspacePool::CreateWorkspace", "GeomVol003",
                  FatalException, "Failed to create workspace.");
    }
    else
    {
      fMyWorkspace= geometryWrk;
    }
  }
  else
  {
    G4Exception("GeometryWorspacePool::CreateWorkspace", "GeomVol003",
                FatalException,
                "Cannot create workspace twice for the same thread.");
    geometryWrk= fMyWorkspace; 
  }
  
  return geometryWrk;
}

// ----------------------------------------------------------------------
// Create it (as above) and use it
//
void G4GeometryWorkspacePool::CreateAndUseWorkspace()
{
  (this->CreateWorkspace())->UseWorkspace();
}

// ----------------------------------------------------------------------
// Reuse an existing workspace - or create a new one if needed.
//
G4GeometryWorkspace* G4GeometryWorkspacePool::FindOrCreateWorkspace()
{
  G4GeometryWorkspace* geometryWrk= fMyWorkspace;
  if( !geometryWrk )
  { 
    geometryWrk= this->CreateWorkspace(); 
  } 
  geometryWrk->UseWorkspace();
  
  fMyWorkspace= geometryWrk; // assign it for use by this thread.
  return geometryWrk;
}

// ----------------------------------------------------------------------
//
void G4GeometryWorkspacePool::Recycle( G4GeometryWorkspace *geometryWrk )
{
  geometryWrk->ReleaseWorkspace(); 
  // if( fWarehouse ){
  // } else {
  delete geometryWrk;
  // }
}

// ----------------------------------------------------------------------
//
void G4GeometryWorkspacePool::CleanUpAndDestroyAllWorkspaces()
{
   if (fMyWorkspace)
   {
      fMyWorkspace->DestroyWorkspace();
      delete fMyWorkspace;
      fMyWorkspace=0;
   }
}

// ----------------------------------------------------------------------
//
G4GeometryWorkspacePool::G4GeometryWorkspacePool()
{
  // fWarehouse=0;
}

// ----------------------------------------------------------------------
//
G4GeometryWorkspacePool::~G4GeometryWorkspacePool()
{
}
