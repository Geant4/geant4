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

//   Create, recycle and Workspaces
//  Only one object of this type exists - it is a true Singleton (not per thread).

#ifndef G4PARTICLESWORKSPACEPOOL_HH
#define G4PARTICLESWORKSPACEPOOL_HH

class G4ParticlesWorkspace;
// #include "G4ParticlesWorkspace.hh"

//Note: Investigate potential to make it a template - on the workspace type

class G4ParticlesWorkspacePool
{
  public: 
     static G4ParticlesWorkspacePool* GetInstance();

   // For use with simple MT mode - each thread gets a workspace and uses it until end
     G4ParticlesWorkspace* CreateWorkspace();

     void CreateAndUseWorkspace();  // Create it (as above) and use it
  
   // For use with 'dynamic' model of threading - workspaces can be recycled
     G4ParticlesWorkspace* FindOrCreateWorkspace();
      // Reuse an existing workspace - or create a new one if needed.
      // This will never fail, except if system is out of resources

     inline G4ParticlesWorkspace* GetWorkspace() { return fMyWorkspace; } 
      // Give back the existing, active workspace for my thread / task

     void Recycle( G4ParticlesWorkspace * );
      // Keep the unused Workspace - for recycling

     void CleanUpAndDestroyAllWorkspaces();
     // To be called once at the end of the job
  
 protected:
      void ReleaseAndDestroyWorkspace(G4ParticlesWorkspace*);
      // Destroy workspace - after releasing it

      // void RegisterWarehouse( G4GeometryWarehouse *) ;
      // The (optional) warehouse keeps a list of free workspaces

 private:
      G4ParticlesWorkspacePool();
     ~G4ParticlesWorkspacePool();

 private: 
     static G4ParticlesWorkspacePool* thePool;

     static G4ThreadLocal G4ParticlesWorkspace*    fMyWorkspace;
     // The thread's workspace - if assigned.
};

#endif //G4GEOMETRYWORKSPACEPOOL_HH
