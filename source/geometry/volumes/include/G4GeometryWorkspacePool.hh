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
// $Id: G4GeometryWorkspacePool.hh 103041 2017-03-10 11:47:01Z gcosmo $
//
//
// ------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
// Create and recycle Workspaces
// Only one object of this type exists; true Singleton (not per thread).

// Author: John Apostolakis (CERN), December 2013
// ------------------------------------------------------------

#ifndef G4GEOMETRYWORKSPACEPOOL_HH
#define G4GEOMETRYWORKSPACEPOOL_HH

#include "G4Types.hh"

class G4GeometryWorkspace;

// Note: Investigate potential to make it a template - on the workspace type

class G4GeometryWorkspacePool
{
  public: 

    static G4GeometryWorkspacePool* GetInstance();

    G4GeometryWorkspace* CreateWorkspace();
      // For use with simple MT mode
      // Each thread gets a workspace and uses it until end

    void CreateAndUseWorkspace();
      // Create it (as above) and use it
  
    G4GeometryWorkspace* FindOrCreateWorkspace();
      // For use with 'dynamic' model of threading. Workspaces can be recycled.
      // Reuse an existing workspace or create a new one if needed.
      // This will never fail, except if system is out of resources

    inline G4GeometryWorkspace* GetWorkspace() { return fMyWorkspace; } 
      // Give back the existing, active workspace for my thread / task

    void Recycle( G4GeometryWorkspace * );
      // Keep the unused Workspace for recycling

    void CleanUpAndDestroyAllWorkspaces();
      // To be called once at the end of the job
  
  protected:

    // void RegisterWarehouse( G4GeometryWarehouse *);
      // The (optional) warehouse keeps a list of free workspaces

  private:

     G4GeometryWorkspacePool();
    ~G4GeometryWorkspacePool();

  private: 

    static G4GeometryWorkspacePool* thePool;

    // void* fWarehouse;
    // G4GeometryWarehouse* fWarehouse;
      // This is "void" to hide the actual container type for workspaces
      // Implementations: a simple STL contaier for MT or something better
      // for tbb.
      // Further development: template parameter with portable default

    static G4ThreadLocal G4GeometryWorkspace* fMyWorkspace;
      // The thread's workspace - if assigned.
      //  --> Can we do without this ?  It is dirty!!
};

#endif
