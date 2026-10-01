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
// G4VoxelisationHelper
//
// Class description:
//
// A helper class to undertake voxelisation in parallel, aiding
// and off-loading the work of G4GeometryManager.
//
// Only one instance must exist, and it must be owned by G4GeometryManager
// -- by the master thread's instance for the time being.
//
// It is designed to be thread safe, and its methods re-entrant.
// The following methods are expected to be called asynchronously
// by threads safely.

// Author: John Apostolakis (CERN), 10.02.2025
// -------------------------------------------------------------------
#ifndef G4VOXELISATIONHELPER_HH
#define G4VOXELISATIONHELPER_HH

#include "G4SmartVoxelStat.hh"
#include "G4Timer.hh"
#include "G4Types.hh"
#include "G4ios.hh"

#include <vector>

class G4VPhysicalVolume;

/**
 * @brief G4VoxelisationHelper is a helper class to undertake voxelisation
 * in parallel, aiding and off-loading the work of G4GeometryManager.
 * @ingroup geometry_management
 */

class G4VoxelisationHelper
{
  public:

    /**
     * Constructor & default destructor.
     */
    G4VoxelisationHelper();
    ~G4VoxelisationHelper() = default;

    /**
     * Key method that creates a list of volumes and resets the state
     * to prepare for the parallel optimisation.
     */
    void PrepareParallelOptimisation(G4bool allOpts, G4bool verbose);

    /**
     * Contributes to voxel optimisation until all work is done.
     * To be called by a worker thread initialisation, not by the user.
     */
    void UndertakeOptimisation();

    // Methods that can be used by this Helper, Geometry Manager and others
    // --------------------------------------------------------------------

    /**
     * Methods to report statistics on the voxelisation process.
     */
    static void ReportVoxelStats(std::vector<G4SmartVoxelStat>& stats, G4double totalCpuTime,
                                 std::ostream& os = G4cout);
    void ReportVoxelInfo(G4LogicalVolume* logVolume, std::ostream& os);

    /**
     * Prints collected statistics once, if verbose optimisation is complete.
     * To be called by the master after synchronising with the workers.
     */
    void ReportParallelOptimisationStats();

    /**
     * Sets voxel statistics verbosity, also retained for subsequent preparations.
     */
    inline void SetVerbosity(G4bool verbose);

    /**
     * Returns true if all workers are finished (or all work is done).
     */
    G4bool IsParallelOptimisationFinished();

    // Auxiliary method - may be useful elsewhere

    /**
     * Check that volumes marked to optimise are done, and report number
     * of those that are missing voxel header.
     */
    G4int CheckOptimisation();

  private:  // Methods used to implement the parallel optimisation

    /**
     * Builds a vector of relevant volumes.
     */
    void CreateListOfVolumesToOptimise(G4bool allOpts, G4bool verbose);

    /**
     * Returns a pointer to the logical volume to optimise.
     */
    G4LogicalVolume* ObtainVolumeToOptimise();

    /**
     * Prepares for doing the work in parallel.
     * Called in preparation (not MT safe).
     */
    void ReSetParallelOptimisation(G4bool verbose);

    /**
     * Resets (empties) the list of candidate volumes for optimisation.
     * Must be called when optimisation is finished.
     */
    void ResetListOfVolumesToOptimise();

    /**
     * Thread-safe method for a worker to report it's finished its work.
     * It counts the number of workers that finished, and returns count.
     * It counts the number of volumes optimised and records completion
     * once all volumes have been optimised.
     */
    G4int ReportWorkerIsDoneOptimising(unsigned int numVolumesOptimised);

    /**
     * Method called once when all volumes are optimised, under the results lock.
     */
    void RecordOptimisationIsFinished(G4bool verbose);

    /**
     * Waits until the voxelisation is all done.
     */
    void WaitForVoxelisationFinish(G4bool verbose = false);

  private:

    /** The list of volumes which threads need to optimise. */
    std::vector<G4LogicalVolume*> fVolumesToOptimise;

    /** Iterator used by UndertakeOptimisation(). */
    std::vector<G4LogicalVolume*>::const_iterator fLogVolumeIterator;

    /** Statistics container shared by all workers. */
    std::vector<G4SmartVoxelStat> fGlobVoxelStats;

    // Flags for parallel initialization
    // ---------------------------------

    /** Explicit verbosity request, independent of each geometry closure. */
    G4bool fRequestedVerboseParallel = false;
    G4bool fVerboseParallel = false;
    G4bool fVoxelStatsReported = false;

    G4bool fParallelVoxelOptimisationUnderway = false;  // It has started
    G4bool fParallelVoxelOptimisationFinished = false;  // It is done

    // Statistics for parallel Optimisation - used in 'verbose' mode
    // ------------------------------------

    /** Counters. */
    G4double fSumVoxelTime = -9999999.9999;
    G4int fNumberThreadsReporting = -99999;  // Must be set correctly later
    unsigned long fTotalNumberVolumesOptimised = -9999999;

    /** Global voxelisation timer, started and stopped under the results lock. */
    G4Timer fWallClockTimer;
};

// Inline methods

inline void G4VoxelisationHelper::SetVerbosity(G4bool verbose)
{
  fRequestedVerboseParallel = verbose;
  fVerboseParallel = verbose;
}

#endif
