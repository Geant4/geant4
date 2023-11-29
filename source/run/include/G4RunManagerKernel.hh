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
// G4RunManagerKernel
//
// Class description:
//
// This is a class for mandatory control of the Geant4 kernel.
//
// This class is constructed by G4RunManager. If a user adopts his/her own
// class instead of G4RunManager, this class must be instantiated by at the
// very beginning of the application and must be deleted at the very end.
// In addition, the following methods must be invoked in the proper order:
//   DefineWorldVolume()
//   InitializePhysics()
//   RunInitialization()
//   RunTermination()
//
// User must provide his/her own classes derived from the following abstract
// class and register it to G4RunManagerKernel:
//   G4VUserPhysicsList - Particle types, Processes and Cuts
//
// G4RunManagerKernel does not have any event loop. Handling of events
// is managed by G4RunManager.

// Author: M.Asai, 1 August 2003
// --------------------------------------------------------------------
#ifndef G4RunManagerKernel_hh
#define G4RunManagerKernel_hh 1

#include "G4EventManager.hh"
#include "globals.hh"

class G4VUserPhysicsList;
class G4VPhysicalVolume;
class G4Region;
class G4ExceptionHandler;
class G4StackManager;
class G4TrackingManager;
class G4PrimaryTransformer;

class G4RunManagerKernel
{
  public:
    // Static method returning the singleton pointer of
    // G4RunManagerKernel or its derived class.
    static G4RunManagerKernel* GetRunManagerKernel();

    // The constructor and the destructor. The user must construct this class
    // object at the beginning of his/her main() and must delete it at the
    // bottom of the main(), unless he/she used G4RunManager.
    G4RunManagerKernel();
    virtual ~G4RunManagerKernel();

    void DefineWorldVolume(G4VPhysicalVolume* worldVol, G4bool topologyIsChanged = true);

    // This method must be invoked if the geometry setup has been changed
    // between runs. The flag "topologyIsChanged" will specify if the
    // geometry topology is different from the original one used in the
    // previous run; if not, it must be set to false, so that the original
    // optimisation and navigation history is preserved. This method is
    // invoked also at initialisation.
    void WorkerDefineWorldVolume(G4VPhysicalVolume* worldVol, G4bool topologyIsChanged = true);

    // This method must be invoked at least once with a valid concrete
    // implementation of user physics list.
    void SetPhysics(G4VUserPhysicsList* uPhys);

    // This method must be invoked at least once to build physics processes.
    void InitializePhysics();

    // Trigger geometry closing and physics table constructions.
    // It returns TRUE if all procedures went well.
    G4bool RunInitialization(G4bool fakeRun = false);

    // Set the application state to 'Idle' so that the user can modify
    // physics/geometry.
    void RunTermination();

    // Update region list. This method is mandatory before invoking the
    // following two dump methods.
    // At RunInitialization(), this method is automatically invoked.
    void UpdateRegion();

    // Dump information of a region.
    void DumpRegion(const G4String& rname) const;

    // Dump information of a region.
    // If the pointer is NULL, all regions are shown.
    void DumpRegion(G4Region* region = nullptr) const;

    void WorkerUpdateWorldVolume();

    // This method must be invoked (or equivalent UI commands can be used)
    // in case the user changes his/her detector geometry.
    // This method is automatically invoked from DefineWorldVolume().
    inline void GeometryHasBeenModified() { geometryNeedsToBeClosed = true; }

    // This method must be invoked in case the user changes his/her physics
    // process(es), e.g. (in)activate some processes. Once this method is
    // invoked, regardless of cuts changed or not, BuildPhysicsTable() of
    // a PhysicsList is invoked for refreshing all physics tables.
    inline void PhysicsHasBeenModified() { physicsNeedsToBeReBuilt = true; }

    inline G4EventManager* GetEventManager() const { return eventManager; }
    inline G4StackManager* GetStackManager() const { return eventManager->GetStackManager(); }
    inline G4TrackingManager* GetTrackingManager() const
    {
      return eventManager->GetTrackingManager();
    }
    inline void SetPrimaryTransformer(G4PrimaryTransformer* pt)
    {
      eventManager->SetPrimaryTransformer(pt);
    }
    inline G4PrimaryTransformer* GetPrimaryTransformer() const
    {
      return eventManager->GetPrimaryTransformer();
    }

    inline const G4String& GetVersionString() const { return versionString; }

    inline void SetVerboseLevel(G4int vl) { verboseLevel = vl; }

    inline void SetGeometryToBeOptimized(G4bool vl)
    {
      if (geometryToBeOptimized != vl) {
        geometryToBeOptimized = vl;
        geometryNeedsToBeClosed = true;
      }
    }

    inline G4int GetNumberOfParallelWorld() const { return numberOfParallelWorld; }
    inline void SetNumberOfParallelWorld(G4int i) { numberOfParallelWorld = i; }

    inline G4VUserPhysicsList* GetPhysicsList() const { return physicsList; }

    inline G4VPhysicalVolume* GetCurrentWorld() const { return currentWorld; }

    inline G4int GetNumberOfStaticAllocators() const { return numberOfStaticAllocators; }

    enum RMKType
    {
      sequentialRMK,
      masterRMK,
      workerRMK
    };

  protected:
    // Constructor to be used by derived classes.
    G4RunManagerKernel(RMKType rmkType);

    // Called by DefineWorldVolume().
    void SetupDefaultRegion();
    void SetupPhysics();
    void ResetNavigator();
    void BuildPhysicsTables(G4bool fakeRun);
    void CheckRegions();

    // This method will setup the G4VProcesses instances to have a reference
    // to the process instance created by the master thread.
    // See G4VProcess::GetMasterProcess().
    virtual void SetupShadowProcess() const;

    void PropagateGenericIonID();

  private:
    void CheckRegularGeometry();
    G4bool ConfirmCoupledTransportation();
    void SetScoreSplitter();

  protected:
    RMKType runManagerKernelType;
    G4Region* defaultRegion = nullptr;
    G4Region* defaultRegionForParallelWorld = nullptr;
    G4bool geometryNeedsToBeClosed = true;

  private:
    G4VUserPhysicsList* physicsList = nullptr;
    G4VPhysicalVolume* currentWorld = nullptr;
    G4bool geometryInitialized = false;
    G4bool physicsInitialized = false;
    G4bool geometryToBeOptimized = true;
    G4bool physicsNeedsToBeReBuilt = true;
    G4int verboseLevel = 0;
    G4int numberOfParallelWorld = 0;

    G4EventManager* eventManager = nullptr;
    G4ExceptionHandler* defaultExceptionHandler = nullptr;
    G4String versionString = "";

    static G4ThreadLocal G4RunManagerKernel* fRunManagerKernel;

    G4int numberOfStaticAllocators = 0;
};

#endif
