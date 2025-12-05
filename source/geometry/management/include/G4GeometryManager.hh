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
// The class is a 'singleton', with access via the static method
// G4GeometryManager::GetInstance().

// Author: Paul Kent (CERN), 26.07.1995 - Initial version
//         John Apostolakis (CERN), 12.06.2024 - Added parallel optimisation
// --------------------------------------------------------------------
#ifndef G4GEOMETRYMANAGER_HH
#define G4GEOMETRYMANAGER_HH

#include <vector>

#include "G4Types.hh"
#include "G4SmartVoxelStat.hh"
#include "G4ios.hh"

class G4VPhysicalVolume;
class G4Timer;
class G4VoxelisationHelper;

/**
 * @brief G4GeometryManager is a singleton class responsible for high level
 * geometrical functions, and for high level objects in the geometry subdomain.
 */

class G4GeometryManager
{
  public:
  
    /**
     * Destructor; called by G4RunManagerKernel.
     */
    ~G4GeometryManager();

    /**
     * Returns a pointer to the singleton instance of the class, creating
     * it if not existing.
     */
    static G4GeometryManager* GetInstance();

    /**
     * Simply returns a pointer to the singleton instance.
     */
    static G4GeometryManager* GetInstanceIfExist();

    /**
     * Closes ('locks') the geometry: performs sanity and 'completion' checks
     * and optionally [default=yes] builds the optimisation structure.
     * Applies to just a specific subtree if a physical volume is specified.
     *  @param[in] pOptimise Flag to enabling/disabling optimisation structure.
     *  @param[in] verbose Flag for verbosity.
     *  @param[in] vol Optional pointer to a physical volume (subtree) for
     *             optimisation.
     *  @returns true if process succeeds.
     */
    G4bool CloseGeometry(G4bool pOptimise = true,
                         G4bool verbose = false,
                         G4VPhysicalVolume* vol = nullptr);
 
    /**
     * Opens ('unlocks') the geometry and removes the optimisation structure if
     * present. Applies to just a specific subtree if a physical volume is
     * specified.
     *  @param[in] vol Optional pointer to a physical volume (subtree) for
     *             optimisation.
     */
    void OpenGeometry(G4VPhysicalVolume* vol = nullptr);

    /**
     * Returns true/false according to the state of the optimised geometry.
     */
    G4bool IsGeometryClosed();

    /**
     * Sets the maximum extent of the world volume.
     * The operation is allowed only if NO solids have been created already.
     */
    void SetWorldMaximumExtent(G4double worldExtent);

    // Methods to choose or undertake Parallel Optimisation
    // ----------------------------------------------------

    /**
     * Requests optimisation using threads (if MT is enabled & used ).
     */
    void OptimiseInParallel(G4bool val = true);
  
    /**
     * Method that contributes to (voxel) optimisation until all work is done.
     * To be called by a worker thread initialisation - not by the user.
     */
    void UndertakeOptimisation();

    /**
     * Detailed method for user to request parallel optimisation (if verbosity
     * is required). Calling this method is enough to ask for it.
     * It will be used if Geant4 is built with MT/tasks.
     * Parallelism will be used if Geant4 is built with MT/tasks.
     */
    void RequestParallelOptimisation(G4bool val = true, G4bool verbose = true);

    /**
     * Checks whether parallel optimisation was requested and configured.
     */
    G4bool IsParallelOptimisationConfigured();

    /**
     * Reports whether parallel optimisation was completed.
     */
    G4bool IsParallelOptimisationFinished();

  private:

    /**
     * Private Constructor. Set the geometry to be open.
     */
    G4GeometryManager();

    // Create the voxel optimisation information. Return whether work is completed
    // (In case of parallel optimisation, it will be done later by the workers)

    /**
     * Optimises all volumes or just multi-volumes (parameterisations, etc. ).
     */
    G4bool BuildOptimisations(G4bool allOpt, G4bool verbose = false);  

    /**
     * Optimises one volume or subtree only.
     */
    G4bool BuildOptimisations(G4bool allOpt, G4VPhysicalVolume* vol);

    /**
     * Builds the optimisation structure in sequential mode.
     */
    void BuildOptimisationsSequential(G4bool allOpts, G4bool verbose = true);

    /**
     * Methods to delete the optimisation structure.
     */
    void DeleteOptimisations();
    void DeleteOptimisations(G4VPhysicalVolume* vol);

  private:

    /** The static instance. */
    static G4ThreadLocal G4GeometryManager* fgInstance;

    /** Pointer to the voxelisation helper. */
    static G4VoxelisationHelper* fParallelVoxeliser;

    G4bool fIsClosed = false;

    // Flags for parallel initialization
    // ---------------------------------
    static G4bool fParallelVoxelOptimisationRequested;
   
    /** Not just requested, but adopted (i.e. also in MT/tasking mode). */
    static G4bool fOptimiseInParallelConfigured;
};

#endif
