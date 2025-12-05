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
// Class G4GeometryManager implementation
//
// Author: Paul Kent (CERN), 26.07.1995 - Initial version
//         John Apostolakis (CERN), 12.06.2024 - Added parallel optimisation
// --------------------------------------------------------------------

#include <iomanip>

#include "G4ios.hh"
#include "G4Timer.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

// Needed for building optimisations
//
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SmartVoxelHeader.hh"
#include "voxeldefs.hh"

// Needed for setting the extent for tolerance value
//
#include "G4GeometryTolerance.hh"
#include "G4SolidStore.hh"
#include "G4VSolid.hh"

// Needed for parallel optimisation
#include "G4AutoLock.hh"
#include "G4VoxelisationHelper.hh"

namespace
{
  // Mutex to obtain a volume to optimise
  G4Mutex createHelperMutex = G4MUTEX_INITIALIZER;
}

// ***************************************************************************
// Static class data
// ***************************************************************************
//
G4ThreadLocal G4GeometryManager* G4GeometryManager::fgInstance = nullptr;

G4VoxelisationHelper* G4GeometryManager::fParallelVoxeliser = nullptr;
 // Only one instance created by the master thread's G4GeometryManager
 // Expected Future: data member (when only one instance)

// Static *global* class data
G4bool G4GeometryManager::fParallelVoxelOptimisationRequested = true;
  // Records User choice to use parallel voxel optimisation (or not)

G4bool G4GeometryManager::fOptimiseInParallelConfigured = false;
  // Configured = requested && available (ie if MT or Threads is used)
  // Value calculated during each effort to optimise

// ***************************************************************************
// Contructor
// ***************************************************************************
//
G4GeometryManager::G4GeometryManager()
{
  fIsClosed = false;

  // Ensure that the helper objects are created
  // even if there is no master thread
  // Note: no longer needed once G4GeometryManager is a canonical Singleton
  G4AutoLock lock(createHelperMutex);
  if( fParallelVoxeliser == nullptr )
  {
    fParallelVoxeliser= new G4VoxelisationHelper();
  }
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4GeometryManager::~G4GeometryManager()
{
  fgInstance = nullptr;
  fIsClosed = false;

  // Only the master thread can delete the helper objects
  //  - a different mechanism is needed for setups with no master
  // TODO: See if shared_ptr could help here?
  if( G4Threading::IsMasterThread() )
  {
    if( fParallelVoxeliser != nullptr )
    {
      delete fParallelVoxeliser;
      fParallelVoxeliser = nullptr;
    }
  }
}

// ***************************************************************************
// Closes geometry - performs sanity checks and optionally builds optimisation
// for placed volumes (always built for replicas & parameterised).
// NOTE: Currently no sanity checks are performed.
// Applies to just a specific subtree if a physical volume is specified.
// ***************************************************************************
//
G4bool G4GeometryManager::CloseGeometry(G4bool pOptimise, G4bool verbose,
                                        G4VPhysicalVolume* pVolume)
{
  if (!fIsClosed && G4Threading::IsMasterThread())
  {
    G4bool workDone= false;
    if (pVolume != nullptr)
    {
      workDone= BuildOptimisations(pOptimise, pVolume);
    }
    else
    {
      workDone= BuildOptimisations(pOptimise, verbose);
    }
    fIsClosed = workDone; // Sequential will be done; parallel ongoing
  }
  return true;
}

// ***************************************************************************
// Inform whether closing of geometry has finished
// ***************************************************************************

G4bool G4GeometryManager::IsGeometryClosed() 
{ 
  if( fOptimiseInParallelConfigured )
  {
    fIsClosed= fParallelVoxeliser->IsParallelOptimisationFinished();
  }
  return fIsClosed;
}

// ***************************************************************************
// Opens the geometry and removes optimisations (optionally, related to just
// the specified logical-volume).
// Applies to just a specific subtree if a physical volume is specified.
// ***************************************************************************
//
void G4GeometryManager::OpenGeometry(G4VPhysicalVolume* pVolume)
{
  if (fIsClosed && G4Threading::IsMasterThread())
  {
    if (pVolume != nullptr)
    {
      DeleteOptimisations(pVolume);
    }
    else
    {
      DeleteOptimisations();
    }
    fIsClosed = false;
    // fGeometryCloseRequested= false;
  }
}

// ***************************************************************************
// Returns the instance of the singleton.
// Creates it in case it's called for the first time.
// ***************************************************************************
//
G4GeometryManager* G4GeometryManager::GetInstance()
{
  if (fgInstance == nullptr)
  {
    fgInstance = new G4GeometryManager;
  }
  return fgInstance;
}

// ***************************************************************************
// Returns the instance of the singleton.
// ***************************************************************************
//
G4GeometryManager* G4GeometryManager::GetInstanceIfExist()
{
  return fgInstance;
}

// ***************************************************************************
// Simplest user method to request parallel optimisation.
// ***************************************************************************
//
void G4GeometryManager::OptimiseInParallel( G4bool val )
{
  RequestParallelOptimisation(val);
}

// ***************************************************************************
// Respond whether parallel optimisation is done
// ***************************************************************************
//
G4bool G4GeometryManager::IsParallelOptimisationFinished()
{ 
  return fParallelVoxeliser->IsParallelOptimisationFinished(); 
}


// ***************************************************************************
// Respond whether parallel optimisation is configured
// ***************************************************************************
// 
G4bool G4GeometryManager::IsParallelOptimisationConfigured()
{
  return fOptimiseInParallelConfigured;
}

// ***************************************************************************
// Creates optimisation info. Builds all voxels if allOpts=true
// otherwise it builds voxels only for replicated volumes.
// ***************************************************************************
// Returns whether optimisation is finished
//
G4bool G4GeometryManager::BuildOptimisations(G4bool allOpts, G4bool verbose)
{
  G4bool finishedOptimisation = false;
  
  fOptimiseInParallelConfigured = fParallelVoxelOptimisationRequested
                               && G4Threading::IsMultithreadedApplication();

  if( fOptimiseInParallelConfigured )
  {
    fParallelVoxeliser->PrepareParallelOptimisation(allOpts, verbose);
  }
  else
  {
    BuildOptimisationsSequential(allOpts, verbose);
    finishedOptimisation= true;
    fIsClosed= true;
  }

  return finishedOptimisation;
}

// ***************************************************************************
// Creates optimisation info. Builds all voxels if allOpts=true
// otherwise it builds voxels only for replicated volumes.
//
// This is the original sequential implementation of this method; was called
// - at first initialisation to create voxels,
// - at re-initialisation if the geometry has changed.
// ***************************************************************************
//
void G4GeometryManager::BuildOptimisationsSequential(G4bool allOpts,
                                                     G4bool verbose)
{
  G4Timer timer;
  G4Timer allTimer;
  std::vector<G4SmartVoxelStat> stats;
  
  if (verbose)  { allTimer.Start(); }
  
  G4LogicalVolumeStore* Store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* volume;
  G4SmartVoxelHeader* head;

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << G4endl
     << "*** G4GeometryManager::BuildOptimisationsSequential() called on tid "
     << G4Threading::G4GetThreadId() << " all-opts= " << allOpts << G4endl;
#endif
  
  for (auto & n : *Store)
  {
    if (verbose) { timer.Start(); }
    volume=n;
    // For safety, check if there are any existing voxels and
    // delete before replacement
    //
    head = volume->GetVoxelHeader();
    delete head;
    volume->SetVoxelHeader(nullptr);
    if (    ( (volume->IsToOptimise())
             && (volume->GetNoDaughters()>=kMinVoxelVolumesLevel1&&allOpts) )
        || ( (volume->GetNoDaughters()==1)
            && (volume->GetDaughter(0)->IsReplicated())
            && (volume->GetDaughter(0)->GetRegularStructureId()!=1) ) )
    {
#ifdef G4GEOMETRY_VOXELDEBUG
    G4cout << "** G4GeometryManager::BuildOptimisationsSequential()"
           << "   Examining logical volume name = '" << volume->GetName()
           << "'  #daughters= " << volume->GetNoDaughters()  << G4endl;
#endif
      head = new G4SmartVoxelHeader(volume);
      
      if (head != nullptr)
      {
        volume->SetVoxelHeader(head);
      }
      else
      {
        std::ostringstream message;
        message << "VoxelHeader allocation error." << G4endl
                << "Allocation of new VoxelHeader" << G4endl
                << "        for volume '" << volume->GetName() << "' failed.";
        G4Exception("G4GeometryManager::BuildOptimisations()", "GeomMgt0003",
                    FatalException, message);
      }
      if (verbose)
      {
        timer.Stop();
        stats.emplace_back( volume, head,
                           timer.GetSystemElapsed(),
                           timer.GetUserElapsed() );
      }
    }
    else
    {
      // Don't create voxels for this node
#ifdef G4GEOMETRY_VOXELDEBUG
      auto numDaughters = volume->GetNoDaughters();
      G4cout << "- Skipping logical volume with " << numDaughters
             << " daughters and name = '" << volume->GetName() << "' " << G4endl;
      if( numDaughters > 1 )
      {
        G4cout << "[Placement]";
      }
      else
      {
        if( numDaughters == 1 )
        {
          G4cout << ( volume->GetDaughter(0)->IsReplicated() ? "[Replicated]"
                                                             : "[Placement]" );
        }
      }
      G4cout << G4endl;
#endif
    }
  }
  if (verbose)
  {
    allTimer.Stop();
    
    G4VoxelisationHelper::ReportVoxelStats( stats, allTimer.GetSystemElapsed()
                                                 + allTimer.GetUserElapsed() );
  }
}

// ***************************************************************************
// Method which user calls to ask for parallel optimisation (or turn it off).
// ***************************************************************************
//
void G4GeometryManager::RequestParallelOptimisation(G4bool flag, G4bool verbose)
{
  fParallelVoxelOptimisationRequested = flag;
  if( flag )
  {
     fParallelVoxeliser->SetVerbosity(verbose);
  }
}

// ***************************************************************************
// Method for a thread/task to contribute dynamically to Optimisation
// ***************************************************************************
//
void G4GeometryManager::UndertakeOptimisation()
{
   fParallelVoxeliser->UndertakeOptimisation();
}


// ***************************************************************************
// Creates Optimisation info for the specified volumes subtree.
// ***************************************************************************
// Returns whether all work is done
//
G4bool G4GeometryManager::BuildOptimisations(G4bool allOpts,
                                             G4VPhysicalVolume* pVolume)
{
  if (pVolume == nullptr) { return false; }

  // Retrieve the mother logical volume, if not NULL,
  // otherwise apply global optimisation for the world volume
  //
  G4LogicalVolume* tVolume = pVolume->GetMotherLogical();
  if (tVolume == nullptr)
  {
    G4bool done=BuildOptimisations(allOpts, false);
    return done;
  }

  G4SmartVoxelHeader* head = tVolume->GetVoxelHeader();
  delete head;
  tVolume->SetVoxelHeader(nullptr);
  if (    ( (tVolume->IsToOptimise())
         && (tVolume->GetNoDaughters()>=kMinVoxelVolumesLevel1&&allOpts) )
       || ( (tVolume->GetNoDaughters()==1)
         && (tVolume->GetDaughter(0)->IsReplicated()) ) ) 
  {
    head = new G4SmartVoxelHeader(tVolume);
    if (head != nullptr)
    {
      tVolume->SetVoxelHeader(head);
    }
    else
    {
      std::ostringstream message;
      message << "VoxelHeader allocation error." << G4endl
              << "Allocation of new VoxelHeader" << G4endl
              << "        for volume " << tVolume->GetName() << " failed.";
      G4Exception("G4GeometryManager::BuildOptimisations()", "GeomMgt0003",
                  FatalException, message);
    }
  }
  else
  {
    // Don't create voxels for this node
#ifdef G4GEOMETRY_VOXELDEBUG
    G4cout << "** G4GeometryManager::BuildOptimisations()" << G4endl
           << "     Skipping logical volume name = " << tVolume->GetName()
           << G4endl;
#endif
  }

  // Scan recursively the associated logical volume tree
  //
  tVolume = pVolume->GetLogicalVolume();
  if (tVolume->GetNoDaughters() != 0)
  {
    BuildOptimisations(allOpts, tVolume->GetDaughter(0));
  }

  return true; // Work is all done -- no parallelism currently
}

// ***************************************************************************
// Removes all optimisation info.
// Loops over all logical volumes, deleting non-null voxels pointers.
// ***************************************************************************
//
void G4GeometryManager::DeleteOptimisations()
{
  G4LogicalVolume* tVolume = nullptr;
  G4LogicalVolumeStore* Store = G4LogicalVolumeStore::GetInstance();
  for (auto & n : *Store)
  {
    tVolume=n;
    delete tVolume->GetVoxelHeader();
    tVolume->SetVoxelHeader(nullptr);
  }
}

// ***************************************************************************
// Removes optimisation info for the specified subtree.
// Scans recursively all daughter volumes, deleting non-null voxels pointers.
// ***************************************************************************
//
void G4GeometryManager::DeleteOptimisations(G4VPhysicalVolume* pVolume)
{
  if (pVolume == nullptr) { return; }

  // Retrieve the mother logical volume, if not NULL,
  // otherwise global deletion to world volume.
  //
  G4LogicalVolume* tVolume = pVolume->GetMotherLogical();
  if (tVolume == nullptr) { return DeleteOptimisations(); }
  delete tVolume->GetVoxelHeader();
  tVolume->SetVoxelHeader(nullptr);

  // Scan recursively the associated logical volume tree
  //
  tVolume = pVolume->GetLogicalVolume();
  if (tVolume->GetNoDaughters() != 0)
  {
    DeleteOptimisations(tVolume->GetDaughter(0));
  }
}

// ***************************************************************************
// Sets the maximum extent of the world volume. The operation is allowed only
// if NO solids have been created already.
// ***************************************************************************
//
void G4GeometryManager::SetWorldMaximumExtent(G4double extent)
{
  if (!G4SolidStore::GetInstance()->empty())
  {
     // Sanity check to assure that extent is fixed BEFORE creating
     // any geometry object (solids in this case)
     //
     G4Exception("G4GeometryManager::SetMaximumExtent()",
                 "GeomMgt0003", FatalException,
                 "Extent can be set only BEFORE creating any geometry object!");
  }
  G4GeometryTolerance::GetInstance()->SetSurfaceTolerance(extent);
}
