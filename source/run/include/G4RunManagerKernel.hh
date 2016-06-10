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
// $Id: G4RunManagerKernel.hh 83384 2014-08-21 14:24:13Z gcosmo $
//
// 

// class description:
//
//     This is a class for mandatory control of GEANT4 kernel. 
//     
//     This class is constructed by G4RunManager. If a user uses his/her own
//     class instead of G4RunManager, this class must be instantiated by
//     him/herself at the very beginning of the application and must be deleted
//     at the very end of the application. Also, following methods must be
//     invoked in the proper order.
//       DefineWorldVolume
//       InitializePhysics
//       RunInitialization
//       RunTermination
// 
//     User must provide his/her own classes derived from the following
//     abstract class and register it to the RunManagerKernel. 
//        G4VUserPhysicsList - Particle types, Processes and Cuts
// 
//     G4RunManagerKernel does not have any eveny loop. Handling of events
//     is managed by G4RunManager.
//

#ifndef G4RunManagerKernel_h
#define G4RunManagerKernel_h 1

class G4VUserPhysicsList;

class G4VPhysicalVolume;
class G4Region;
class G4ExceptionHandler;
class G4StackManager;
class G4TrackingManager;
class G4PrimaryTransformer;

#include "globals.hh"
#include "G4EventManager.hh"

class G4RunManagerKernel
{
  public: // with description
    static G4RunManagerKernel* GetRunManagerKernel();
    //  Static method which returns the singleton pointer of G4RunManagerKernel or
    // its derived class.

  private:
    static G4ThreadLocal G4RunManagerKernel* fRunManagerKernel;

  public: // with description

    G4RunManagerKernel();
    virtual ~G4RunManagerKernel();
    //  The constructor and the destructor. The user must construct this class
    // object at the beginning of his/her main() and must delete it at the 
    // bottom of the main(), unless he/she used G4RunManager.
  public:
    enum RMKType { sequentialRMK, masterRMK, workerRMK };
  protected:
    //Constructor to be used by derived classes 
    G4RunManagerKernel(RMKType rmkType);
    RMKType runManagerKernelType;

  public: // with description
    void DefineWorldVolume(G4VPhysicalVolume * worldVol,
                           G4bool topologyIsChanged=true);

    void WorkerDefineWorldVolume(G4VPhysicalVolume * worldVol,
                                G4bool topologyIsChanged=true);

    //  This method must be invoked if the geometry setup has been changed between
    // runs. The flag 'topologyIsChanged' will specify if the geometry topology is
    // different from the original one used in the previous run; if not, it must be
    // set to false, so that the original optimisation and navigation history is
    // preserved. This method is invoked also at initialisation.

    void SetPhysics(G4VUserPhysicsList* uPhys);
    //  This method must be invoked at least once by the user with a valid
    // concrete implementation of user physics list. 

    void InitializePhysics();
    //  This method must be invoked at least once by the user to build physics
    // processes.

    G4bool RunInitialization(G4bool fakeRun=false);
    //  Trigger geometry closing and physics table constructions.
    // It returns TRUE if all procedures went well.

    void RunTermination();
    //  Set the application state to G4State_Idle so that the user can modify
    // physics/geometry.

  public:
    void WorkerUpdateWorldVolume();

  protected:
    void SetupDefaultRegion();
    //Called by DefineWorldVolume
    void SetupPhysics();
    void ResetNavigator();
    void BuildPhysicsTables(G4bool fakeRun);
    void CheckRegions();

  public: // with description
    void UpdateRegion();
    // Update region list. 
    // This method is mandatory before invoking following two dump methods.
    // At RunInitialization(), this method is automatically invoked, and thus
    // the user needs not invoke.

    void DumpRegion(const G4String& rname) const;
    // Dump information of a region.

    void DumpRegion(G4Region* region=0) const;
    // Dump information of a region.
    // If the pointer is NULL, all regions are shown.

  private:
    G4VUserPhysicsList * physicsList;
    G4VPhysicalVolume* currentWorld;
    G4bool geometryInitialized;
    G4bool physicsInitialized;
    G4bool geometryToBeOptimized;
    G4bool physicsNeedsToBeReBuilt;
    G4int verboseLevel;
    G4int numberOfParallelWorld;

    G4EventManager * eventManager;
    G4ExceptionHandler* defaultExceptionHandler;
    G4String versionString;
  protected:
    G4Region* defaultRegion;
    G4Region* defaultRegionForParallelWorld;
    G4bool geometryNeedsToBeClosed;
  public: // with description
    inline void GeometryHasBeenModified()
    { geometryNeedsToBeClosed = true; }
    //  This method must be invoked (or equivalent UI commands can be used)
    // in case the user changes his/her detector geometry.
    // This method is automatically invoked from DefineWorldVolume() method.

    inline void PhysicsHasBeenModified()
    { physicsNeedsToBeReBuilt = true; }
    //  This method must be invoked in case the user changes his/her physics
    // process(es), e.g. (in)activate some processes. Once this method is
    // invoked, regardless of cuts are changed or not, BuildPhysicsTable()
    // of PhysicsList is invoked for refreshing all physics tables.

  public:
    inline G4EventManager* GetEventManager() const
    { return eventManager; }
    inline G4StackManager* GetStackManager() const
    { return eventManager->GetStackManager(); }
    inline G4TrackingManager* GetTrackingManager() const
    { return eventManager->GetTrackingManager(); }
    inline void SetPrimaryTransformer(G4PrimaryTransformer* pt)
    { eventManager->SetPrimaryTransformer(pt); }
    inline G4PrimaryTransformer* GetPrimaryTransformer() const
    { return eventManager->GetPrimaryTransformer(); }

    inline const G4String& GetVersionString() const
    { return versionString; }

    inline void SetVerboseLevel(G4int vl)
    { verboseLevel = vl; }

    inline void SetGeometryToBeOptimized(G4bool vl)
    { 
      if(geometryToBeOptimized != vl)
      {
        geometryToBeOptimized = vl;
        geometryNeedsToBeClosed = true;
      }
    }

    inline G4int GetNumberOfParallelWorld() const
    { return numberOfParallelWorld; }
    inline void SetNumberOfParallelWorld(G4int i)
    { numberOfParallelWorld = i; }

    inline G4VUserPhysicsList* GetPhysicsList() const
    { return physicsList; }

    inline G4VPhysicalVolume* GetCurrentWorld() const
    { return currentWorld; }
  private:
    void CheckRegularGeometry();
    G4bool ConfirmCoupledTransportation();
    void SetScoreSplitter();

    G4int numberOfStaticAllocators;

  public:
    inline G4int GetNumberOfStaticAllocators() const
    { return numberOfStaticAllocators; }
protected:
    virtual void SetupShadowProcess() const;
    // This method will setup the G4VProcesses
    // instances to have a reference to the process instance
    // created by the master thread. See G4VProcess::GetMasterProcess

    void PropagateGenericIonID();
};

#endif

