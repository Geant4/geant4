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
// G4GeometryManager
//
// Class description:
//
// A class responsible for high level geometrical functions, and for
// high level objects in the geometry subdomain.
// The class is a `singleton', with access via the static method
// G4GeometryManager::GetInstance().
//
// Member data:
//
//   - fgInstance
//     Ptr to the unique instance of class (per Thread)

// 26.07.95, P.Kent - Initial version, including optimisation build
// 12.06.24, J.Apostolakis - Added parallel optimisation in workers
// --------------------------------------------------------------------
#ifndef G4GEOMETRYMANAGER_HH
#define G4GEOMETRYMANAGER_HH 1

#include <vector>

#include "G4Types.hh"
#include "G4SmartVoxelStat.hh"
#include "G4ios.hh"

class G4VPhysicalVolume;
class G4Timer;

class G4GeometryManager
{
  public:
  
    G4bool CloseGeometry(G4bool pOptimise = true, G4bool verbose = false,
                         G4VPhysicalVolume* vol = nullptr);
      // Close (`lock') the geometry: perform sanity and `completion' checks
      // and optionally [default=yes] build optimisation information.
      // Applies to just a specific subtree if a physical volume is specified.

    void OpenGeometry(G4VPhysicalVolume* vol = nullptr);
      // Open (`unlock') the geometry and remove optimisation information if
      // present. Applies to just a specific subtree if a physical volume is
      // specified.

    static G4bool IsGeometryClosed();
      // Return true/false according to state of optimised geoemtry.

    void SetWorldMaximumExtent(G4double worldExtent);
      // Set the maximum extent of the world volume. The operation is
      // allowed only if NO solids have been created already.

    static G4GeometryManager* GetInstance();
      // Return ptr to singleton instance of the class, creating it if
      // not existing.

    static G4GeometryManager* GetInstanceIfExist();
      // Return ptr to singleton instance.

    static void OptimizeInParallel(G4bool val = true);
      // Request optimization using threads (if MT is enabled & used ).
  
    void UndertakeOptimisation();
      // Method that contributes to (Voxel) optimisation until all work is done.
      // Must be called by Worker thread initialisation - not a user callable
      // method.

    static void RequestParallelOptimisation(G4bool val = true,
                                            G4bool verbose = true);
      // Detailed method for user to request parallel Optimisation
      // (if verbosity is required). Calling this is enough to ask for it.
      // It will be used if Geant4 is built with MT/tasks.

    static void ChooseSequentialOptimisation(G4bool verbose = false);
      // Simple way to avoid parallel optimisation.
  
    static G4bool IsParallelOptimisationConfigured();
      // Check whether parallel optimisation was requested.
    static G4bool IsParallelOptimisationFinished();
      // Report whether parallel optimisation is done.
  
   ~G4GeometryManager();
      // Destructor; called by G4RunManagerKernel.

  private:

    G4GeometryManager() = default;
      // Private constructor. Set the geometry to be open.

    G4bool BuildOptimisations(G4bool allOpt, G4bool verbose = false);  
       // Optimise all or just multi-volumes (parameterisations, .. ).
    void BuildOptimisations(G4bool allOpt, G4VPhysicalVolume* vol);
       // Optimise one volume or subtree only.
    void DeleteOptimisations();
    void DeleteOptimisations(G4VPhysicalVolume* vol);
  
    void ReportVoxelStats( std::vector<G4SmartVoxelStat>& stats,
                           G4double totalCpuTime,
                           std::ostream &os = G4cout );
    void ReportVoxelInfo(G4LogicalVolume * logVolume, std::ostream& os);
   
    void PrepareParallelOptimisation(G4bool allOpts, G4bool verbose = true);
    void BuildOptimisationsSequential(G4bool allOpts, G4bool verbose = true);

    // Methods for parallel initialization
    void CreateListOfVolumesToOptimise(G4bool allOpts, G4bool verbose);
      // Build vector of relevant volumes.
    G4LogicalVolume* ObtainVolumeToOptimize();

    static G4ThreadLocal G4GeometryManager* fgInstance;
    static G4ThreadLocal G4bool fIsClosed;

    static std::vector<G4LogicalVolume*> fVolumesToOptimize;
      // The list of volumes which threads need to optimize.
    static std::vector<G4LogicalVolume*>::iterator fLogVolumeIterator;
      // Iterator used by UndertakeOptimisation().

    static std::vector<G4SmartVoxelStat> fGlobVoxelStats;
      // Statistics container shared by all workers
  
    static void ConfigureParallelOptimisation(G4bool verbose);
      // Prepare for parallel optimisation.

    G4int ReportWorkerIsDoneOptimising(unsigned int numVolumesOptimized);
      // Thread-safe method for worker to report it's finished its work.
      // It counts the number of workers that finished, and returns count.
      // It counts the number of volumes optimised; if all workers have
      // reported, it results in a 'Finished' state.
  
    static void InformOptimisationIsFinished(G4bool verbose);
      // Returns true if all workers are finished (or all work is done).
  
    static void  ResetListOfVolumesToOptimise();
      // Resets (empties) the list of candidate volumes for optimisation.
      // Must be called when Optimisation is finished.
  
    G4int CheckOptimisation();
      // Check volumes marked to optimised are done, and report number
      // that are missing voxel header.
  
    void WaitForVoxelisationFinish(G4bool verbose = false);
      // Wait until the voxelisation is all done.
  
  private:

    // Flags for parallel initialization
    // ---------------------------------
    static G4bool fVerboseParallel;
    static G4bool fParallelVoxelOptimisationRequested;
      // Flag to register it was requested.
    static G4bool fOptimizeInParallelConfigured;
      // Not just requested, but adopted (i.e. also in MT/tasking mode).
    static G4bool fParallelVoxelOptimisationUnderway; // It has started
    static G4bool fParallelVoxelOptimisationFinished; // It is done
    static G4bool fUsingExistingWorkers;
      // Fact: can and will use existing MT/tasks.

    // Statistics for parallel Optimisation - used in 'verbose' mode
    // ------------------------------------
    static G4double fSumVoxelTime;
    static G4int fNumberThreadsReporting;
    static unsigned int fTotalNumberVolumesOptimized;
      // Counters.
  
    // For Wall Clock time in parallel mode ...
    //
    static G4Timer* fWallClockTimer;   // Owned by master thread
    static G4bool fWallClockStarted;
};

#endif
