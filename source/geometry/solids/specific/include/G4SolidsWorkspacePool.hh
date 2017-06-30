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

#ifndef G4SOLIDSWORKSPACEPOOL_HH
#define G4SOLIDSWORKSPACEPOOL_HH

class G4SolidsWorkspace;

class G4SolidsWorkspacePool
{
  public: 
     static G4SolidsWorkspacePool* GetInstance();

     G4SolidsWorkspace* CreateWorkspace();
       // For use with simple MT mode - each thread gets a workspace
       // and uses it until end

     void CreateAndUseWorkspace();  // Create it (as above) and use it
  
     G4SolidsWorkspace* FindOrCreateWorkspace();
      // For use with 'dynamic' model of threading - workspaces can be recycled
      // Reuse an existing workspace - or create a new one if needed.
      // This will never fail, except if system is out of resources

     inline G4SolidsWorkspace* GetWorkspace() { return fMyWorkspace; } 
      // Give back the existing, active workspace for my thread / task

     void Recycle( G4SolidsWorkspace * );
      // Keep the unused Workspace - for recycling

     void CleanUpAndDestroyAllWorkspaces();
     // To be called once at the end of the job
  
 protected:

      // void RegisterWarehouse( G4GeometryWarehouse *) ;
      // The (optional) warehouse keeps a list of free workspaces

 private:
      G4SolidsWorkspacePool();
     ~G4SolidsWorkspacePool();

 private: 
     static G4SolidsWorkspacePool* thePool;

     static G4ThreadLocal G4SolidsWorkspace* fMyWorkspace;
     // The thread's workspace - if assigned.
};

#endif //G4GEOMETRYWORKSPACEPOOL_HH
