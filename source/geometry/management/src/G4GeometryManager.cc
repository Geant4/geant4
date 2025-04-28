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
// 26.07.95, P.Kent - Initial version, including optimisation build
// 12.06.24, J.Apostolakis - Added parallel optimisation in workers
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

namespace  // Data structures / mutexes for parallel optimisation
{
  // Mutex to obtain a volume to optimise
  G4Mutex obtainVolumeMutex = G4MUTEX_INITIALIZER;

  // Mutex to lock saving of voxel statistics
  G4Mutex voxelStatsMutex = G4MUTEX_INITIALIZER;

  // Mutex to provide Statistics Results
  G4Mutex statResultsMutex = G4MUTEX_INITIALIZER;

  // Mutex to start wall clock (global) timer
  G4Mutex wallClockTimerMutex = G4MUTEX_INITIALIZER;

  // Mutex to write debug output
  G4Mutex outputDbgMutex = G4MUTEX_INITIALIZER;
}

// ***************************************************************************
// Static class data
// ***************************************************************************
//
G4ThreadLocal G4GeometryManager* G4GeometryManager::fgInstance = nullptr;

// Static *global* class data
G4bool G4GeometryManager::fParallelVoxelOptimisationRequested = false;
  // Records User choice to use parallel voxel optimisation (or not)

G4bool  G4GeometryManager::fOptimiseInParallelConfigured = false;
  // Configured = requested && available (ie if MT or Threads is used)
  // Value calculated during each effort to optimise

std::vector<G4LogicalVolume*> G4GeometryManager::fVolumesToOptimise;
std::vector<G4LogicalVolume*>::const_iterator G4GeometryManager::fLogVolumeIterator;

std::vector<G4SmartVoxelStat> G4GeometryManager::fGlobVoxelStats;
// Container for statistics

// Derived state
G4bool G4GeometryManager::fVerboseParallel = false;
G4bool G4GeometryManager::fParallelVoxelOptimisationUnderway = false;
G4bool G4GeometryManager::fParallelVoxelOptimisationFinished = false;
G4double G4GeometryManager::fSumVoxelTime = 0.0;

G4int G4GeometryManager::fNumberThreadsReporting = 0;
unsigned int G4GeometryManager::fTotalNumberVolumesOptimised = 0U;

// For Wall clock
G4Timer* G4GeometryManager::fWallClockTimer = nullptr;
G4bool G4GeometryManager::fWallClockStarted = false;

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4GeometryManager::~G4GeometryManager()
{
  fgInstance = nullptr;
  fIsClosed = false;
  
  if( (fWallClockTimer != nullptr) && G4Threading::IsMasterThread() )
  {
    delete fWallClockTimer;
    fWallClockTimer= nullptr;
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
    if (pVolume != nullptr)
    {
      BuildOptimisations(pOptimise, pVolume);
    }
    else
    {
      BuildOptimisations(pOptimise, verbose);
    }
    fIsClosed = true;
  }
  return true;
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
    
    if( (fWallClockTimer == nullptr) && G4Threading::IsMasterThread() )
    {
      fWallClockTimer = new G4Timer;
    }
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
// Report about Voxel(isation) of a logical volume.
// ***************************************************************************
//
void
G4GeometryManager::ReportVoxelInfo(G4LogicalVolume* logVolume, std::ostream& os)
{
  G4SmartVoxelHeader* head = logVolume->GetVoxelHeader();
  if( head != nullptr )
  {
    os << "** Created optimisations for logical-volume '" 
       << std::setw(50) << logVolume->GetName() << "'" << G4endl
       << "- Result VoxelInfo - START: " << " ptr= " << head << G4endl
       << *head
       << "- Result VoxelInfo -   END. " << G4endl;
  }
  else
  {
    os << "** No optimisation for log-vol " << logVolume->GetName() << G4endl;
  }
  os << "*** Report Voxel Info: END " <<  G4endl;
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
    PrepareParallelOptimisation(allOpts, verbose);
  }
  else
  {
    BuildOptimisationsSequential(allOpts, verbose);
    finishedOptimisation= true;
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
    
    ReportVoxelStats( stats, allTimer.GetSystemElapsed()
                     + allTimer.GetUserElapsed() );
  }
}

// ***************************************************************************
// Creates a list of logical volumes which will be optimised
//    if allOpts=true it lists all voxels
//    otherwise       it lists only the voxels of replicated volumes.
// This list will be used subsequently to build their voxels.
//
// Note: this method is NOT thread safe!
//    It expects to be called only once in each (re)initalisation
//    i.e. either by master thread or a selected thread.
// ***************************************************************************
//
void
G4GeometryManager::CreateListOfVolumesToOptimise(G4bool allOpts, G4bool verbose)
{
  // Prepare the work - must be called only in one thread !!

  G4LogicalVolumeStore* Store = G4LogicalVolumeStore::GetInstance();

  if( !fVolumesToOptimise.empty() )
  {
    ResetListOfVolumesToOptimise();
  }
  
  for (auto & n : *Store)
  {
    G4LogicalVolume* volume=n;
    
    if (    ( (volume->IsToOptimise())
             && (volume->GetNoDaughters()>=kMinVoxelVolumesLevel1&&allOpts) )
        || ( (volume->GetNoDaughters()==1)
            && (volume->GetDaughter(0)->IsReplicated())
            && (volume->GetDaughter(0)->GetRegularStructureId()!=1) ) )
    {
      fVolumesToOptimise.push_back(volume);

      // For safety, must check (later) if there are any existing voxels and
      //     delete before replacement:
      // All 'clients' of this code must do the following:
      //   delete volume->GetVoxelHeader();
      //   volume->SetVoxelHeader(nullptr);
      
#ifdef G4GEOMETRY_VOXELDEBUG
      G4cout << "- Booking  logical volume with " << volume->GetNoDaughters()
      << " daughters and name = '" << volume->GetName() << "' "
      << " -- for optimisation (ie voxels will be built for it). " << G4endl;
#endif
    }
    else
    {
#ifdef G4GEOMETRY_VOXELDEBUG
      G4cout << "- Skipping logical volume with " << volume->GetNoDaughters()
      << " daughters and name = '" << volume->GetName() << "' " << G4endl;
#endif
    }
  }
  
  if(verbose)
  {
    G4cout << "** G4GeometryManager::CreateListOfVolumesToOptimise: "
           << "  Number of volumes for voxelisation = "
           << fVolumesToOptimise.size() << G4endl;
  }
  
  fLogVolumeIterator = fVolumesToOptimise.cbegin();
}

// ***************************************************************************
// Obtain a logical volume from the list of volumes to optimise
// Must be thread-safe: its role is to be called in parallel by threads/tasks!
// Critical method for parallel optimisation - must be correct and fast.
// ***************************************************************************
//
G4LogicalVolume* G4GeometryManager::ObtainVolumeToOptimise()
{
  G4LogicalVolume* logVolume = nullptr;

  G4AutoLock lock(obtainVolumeMutex);
    
  if( fLogVolumeIterator != fVolumesToOptimise.cend() )
  {
    logVolume = *fLogVolumeIterator;
    ++fLogVolumeIterator;
  }
  return logVolume;
}

// ***************************************************************************
// Thread-safe method to clear the list of volumes to Optimise.
// ***************************************************************************
//
void G4GeometryManager::ResetListOfVolumesToOptimise()
{
  G4AutoLock lock(obtainVolumeMutex);

  std::vector<G4LogicalVolume*>().swap(fVolumesToOptimise);
  // Swapping with an empty vector in order to empty it
  // without calling destructors of logical volumes.
  // Must not call clear: i.e. fVolumesToOptimise.clear();

  assert(fVolumesToOptimise.empty());
  fLogVolumeIterator = fVolumesToOptimise.cbegin();
  
  fGlobVoxelStats.clear();
  // Reset also the statistics of volumes -- to avoid double recording.
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
    ConfigureParallelOptimisation(verbose);
  }
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
// Setup up state to enable parallel optimisation by workers.
// ***************************************************************************
//
void G4GeometryManager::ConfigureParallelOptimisation(G4bool verbose)
{
  static G4ThreadLocal unsigned int numCallsConfig = 0;
  ++numCallsConfig;

  if(verbose)
  {
    G4cout << "** G4GeometryManager::ConfigureParallelOptimisation() called "
    << " for the " << numCallsConfig << " time. "
    << " LEAVING all the work (of voxel optimisation) to the threads/tasks !"
    << G4endl;
  }
  fParallelVoxelOptimisationRequested = true;
  fParallelVoxelOptimisationUnderway = false;
  fParallelVoxelOptimisationFinished = false;
  
  // Keep values of options / verbosity for use in threads
  fVerboseParallel = verbose;
  
  // New effort -- reset the total time -- and number of threads reporting
  fSumVoxelTime = 0.0;
  fNumberThreadsReporting = 0;
  fTotalNumberVolumesOptimised = 0;   // Number of volumes done
  
  fWallClockStarted = false;  // Will need to restart it!

  fLogVolumeIterator = fVolumesToOptimise.cbegin();
   // Reset the iterator -- else to be sure that work is done correctly
}

// ***************************************************************************
// Build voxel optimisation in parallel -- prepare the work for threads/tasks
// ***************************************************************************
//
void
G4GeometryManager::PrepareParallelOptimisation(G4bool allOpts, G4bool verbose)
{
  if( verbose )
  {
    G4cout << "** G4GeometryManager::PrepareParallelOptimisation() called."
           << G4endl;
  }
  CreateListOfVolumesToOptimise(allOpts, verbose);
  ConfigureParallelOptimisation(verbose);
}

// ***************************************************************************
// Method for a thread/task to contribute dynamically to Optimisation
// ***************************************************************************
//
void G4GeometryManager::UndertakeOptimisation()
{
  G4bool verbose = fVerboseParallel;
  G4LogicalVolume* logVolume = nullptr;

  fParallelVoxelOptimisationUnderway = true;

  // Start timer - if not already done
  if( ( !fWallClockStarted ) && verbose )
  {
    G4AutoLock startTimeLock(wallClockTimerMutex);
    if( !fWallClockStarted )
    {
      fWallClockTimer->Start();
      fWallClockStarted= true;
    }
  }

  G4Timer fetimer;
  unsigned int numVolumesOptimised = 0;
  
  while( (logVolume = ObtainVolumeToOptimise()) != nullptr )
  {
    if (verbose) { fetimer.Start(); }

    G4SmartVoxelHeader* head = logVolume->GetVoxelHeader();
    delete head;
    logVolume->SetVoxelHeader(nullptr);

    head = new G4SmartVoxelHeader(logVolume);
    //     *********************************
    logVolume->SetVoxelHeader(head);
    
    if (head != nullptr)
    {
      ++numVolumesOptimised;
    }
    else
    {
      G4ExceptionDescription message;
      message << "VoxelHeader allocation error." << G4endl
              << "Allocation of new VoxelHeader" << G4endl
              << "for logical volume " << logVolume->GetName() << " failed.";
      G4Exception("G4GeometryManager::BuildOptimisationsParallel()",
                  "GeomMgt0003", FatalException, message);
    }

    if(verbose)
    {
      fetimer.Stop();
      auto feRealElapsed = fetimer.GetRealElapsed();
      // Must use 'real' elapsed time -- cannot trust user/system time
      // (it accounts for all threads)
      
      G4AutoLock lock(voxelStatsMutex);
      fGlobVoxelStats.emplace_back( logVolume, head,
                          0.0,             // Cannot estimate system time
                          feRealElapsed ); // Use real time instead of user time
      fSumVoxelTime += feRealElapsed;
    }
  }

  G4bool allDone = false;
  G4int myCount= -1;

  myCount = ReportWorkerIsDoneOptimising(numVolumesOptimised);
  allDone = IsParallelOptimisationFinished();

  if( allDone && (myCount == G4Threading::GetNumberOfRunningWorkerThreads()) )
  {
    G4int badVolumes = CheckOptimisation(); // Check all voxels are created!
    if( badVolumes > 0 )
    {
      G4ExceptionDescription errmsg;
      errmsg <<" Expected that all voxelisation work is done, "
             << "but found that voxels headers are missing in "
             << badVolumes << " volumes.";
      G4Exception("G4GeometryManager::UndertakeOptimisation()",
                  "GeomMgt002", FatalException, errmsg);
    }
    
    // Create report

    if( verbose )
    {
      fWallClockTimer->Stop();

      std::ostream& report_stream = std::cout; // G4cout; does not work!
      report_stream << G4endl
        << "G4GeometryManager::UndertakeOptimisation"
        << " -- Timing for Voxel Optimisation" << G4endl;
      report_stream << "  - Elapsed time (real) = " << std::setprecision(4)
        << fWallClockTimer->GetRealElapsed() << "s (wall clock)"
        << ", user " << fWallClockTimer->GetUserElapsed() << "s"
        << ", system " << fWallClockTimer->GetSystemElapsed() << "s."
        << G4endl;
      report_stream << "  - Sum voxel time (real) = " << fSumVoxelTime
                    << "s.";
      report_stream << std::setprecision(6) << G4endl << G4endl;

      ReportVoxelStats( fGlobVoxelStats, fSumVoxelTime, report_stream );
      report_stream.flush();
    }
  }
  else
  {
    WaitForVoxelisationFinish(false);
  }
}

// ***************************************************************************
// Ensure that all the work of voxelisation is done.
// Can be called in GeometryManager methods or externally.
// ***************************************************************************
//
void G4GeometryManager::WaitForVoxelisationFinish(G4bool verbose)
{
  // Must wait until all workers are done ...
  using namespace std::chrono_literals;
  unsigned int trials = 0;
  auto tid = G4Threading::G4GetThreadId();
  
  std::ostream& out_stream = std::cout; // G4cout; does not work!
  while( ! IsParallelOptimisationFinished() )
  {
    // Each thread must wait until all are done ...
    std::this_thread::sleep_for(250ms);
    ++trials;
  }
  
  if( verbose )
  {
    G4AutoLock lock(outputDbgMutex);
    out_stream << G4endl
               << "** UndertakeOptimisation done on tid= " << tid
               <<  " after waiting for " << trials << " trials." << G4endl;
    out_stream.flush();
  }
}

// ***************************************************************************
// Ensure that all logical volumes in list have a voxel-header.
// ***************************************************************************
//
G4int G4GeometryManager::CheckOptimisation()
{
  unsigned int numErrors = 0;
  for ( const auto& logical : fVolumesToOptimise )
  {
    if( logical->GetVoxelHeader() == nullptr ) { ++numErrors; }
  }
  return numErrors;
}

// ***************************************************************************
// Report that current thread/task is done optimising.
// A thread call this method to reports that is is done (finished), and how
// many volumes it optimised. The method:
//   - increments the count of workers that have finished, and return it;
//   - keeps count of number of volumes optimised;
//   - if all works is done (ie all workers have reported) it will result
//     in the 'Finished' state.
// ***************************************************************************
//
G4int
G4GeometryManager::ReportWorkerIsDoneOptimising(unsigned int nVolOptimised)
{
  // Check that all are done and, if so, signal that optimisation is finished
  G4int orderReporting;
  
  G4AutoLock lock(statResultsMutex);
  orderReporting = ++fNumberThreadsReporting;
  fTotalNumberVolumesOptimised += nVolOptimised;
  
  if( fVerboseParallel )
  {
    std::cout << "G4GeometryManager: the " << orderReporting
              << " worker has finished.  "
              <<  "  Total volumes voxelised = "
              << fTotalNumberVolumesOptimised
              <<  " out of " << fVolumesToOptimise.size() << G4endl;
  }

  if ( fNumberThreadsReporting == G4Threading::GetNumberOfRunningWorkerThreads()
       || fTotalNumberVolumesOptimised == fVolumesToOptimise.size() )
  {
    const auto TotalThreads = G4Threading::GetNumberOfRunningWorkerThreads();
    auto tid= G4Threading::G4GetThreadId();

    // -- Some Checks
    if( fTotalNumberVolumesOptimised != fVolumesToOptimise.size() )
    {
      G4ExceptionDescription errmsg;
      errmsg << " [thread " << tid << " ] "
             << " ERROR: Number of volumes 'voxelised' = "
             << fTotalNumberVolumesOptimised
             << " is not equal to the total number requested "
             << fVolumesToOptimise.size() << " !! " << G4endl;
      G4Exception("G4GeometryManager::ReportWorkerIsDoneOptimising()",
                  "GeomMgt0003", FatalException, errmsg);
    }
    
    if( fVerboseParallel )
    {
      if( fNumberThreadsReporting > TotalThreads )
      {
        G4ExceptionDescription errmsg;
        errmsg << " [thread " << tid << " ] "
               << " WARNING: Number of threads 'reporting' = "
               << fNumberThreadsReporting
               << " exceeds the total number of threads "
               << TotalThreads << " !! " << G4endl
               << " *Missed* calling ConfigureParallelOptimisation() to reset. ";
        G4Exception("G4GeometryManager::ReportWorkerIsDoneOptimising()",
                    "GeomMgt1002", JustWarning, errmsg);
      }
      else
      {
        if( fTotalNumberVolumesOptimised == fVolumesToOptimise.size()
            && ( fNumberThreadsReporting < TotalThreads ) )
        {
           G4ExceptionDescription errmsg;
           errmsg << " [thread " << tid << " ] "
                  << " WARNING: All volumes optimised, yet only "
                  << fNumberThreadsReporting << " threads reported out of "
                  << TotalThreads;
           G4Exception("G4GeometryManager::ReportWorkerIsDoneOptimising()",
                       "GeomMgt1002", JustWarning, errmsg);
        }
      }
    }
    // -- End of Checks

    // Report the end
    if( fNumberThreadsReporting <= G4Threading::GetNumberOfRunningWorkerThreads() )
    {
      // Close the work, and (if verbosity is on) report statistics
      InformOptimisationIsFinished(fVerboseParallel);
    }
  }
  
  return orderReporting;
}

// ***************************************************************************
// Inform that all work for parallel optimisation is finished.
// ***************************************************************************
//
void G4GeometryManager::InformOptimisationIsFinished(G4bool verbose)
{
  if(verbose)   // G4cout does not work!
  {
    std::cout << "** G4GeometryManager: All voxel optimisation work is completed!"
              << G4endl;
    std::cout << "   Total number of volumes optimised = "
              << fTotalNumberVolumesOptimised 
              << " of " << fVolumesToOptimise.size() << " expected\n";
    std::cout << "   Number of workers reporting       = "
              << fNumberThreadsReporting
              << " of " << G4Threading::GetNumberOfRunningWorkerThreads()
              << " expected\n";
  }
  assert ( fTotalNumberVolumesOptimised == fVolumesToOptimise.size() );

  fParallelVoxelOptimisationFinished  = true;
  // fParallelVoxelOptimisationRequested = false; // Maintain request for next one!
  fParallelVoxelOptimisationUnderway  = false; // It's no longer underway!
}

// ***************************************************************************
// Creates Optimisation info for the specified volumes subtree.
// ***************************************************************************
//
void G4GeometryManager::BuildOptimisations(G4bool allOpts,
                                           G4VPhysicalVolume* pVolume)
{
  if (pVolume == nullptr) { return; }

  // Retrieve the mother logical volume, if not NULL,
  // otherwise apply global optimisation for the world volume
  //
  G4LogicalVolume* tVolume = pVolume->GetMotherLogical();
  if (tVolume == nullptr)
  {
    BuildOptimisations(allOpts, false);
    return;
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

// ***************************************************************************
// Reports statistics on voxel optimisation when closing geometry.
// ***************************************************************************
//
void
G4GeometryManager::ReportVoxelStats( std::vector<G4SmartVoxelStat> & stats,
                                     G4double totalCpuTime, std::ostream &os )
{
  os << "--------------------------------------------------------------------------------"
     << G4endl;
  os << "G4GeometryManager::ReportVoxelStats -- Voxel Statistics"
         << G4endl << G4endl;
 
  //
  // Get total memory use
  //
  G4int i, nStat = (G4int)stats.size();
  G4long totalMemory = 0;
 
  for( i=0; i<nStat; ++i )  { totalMemory += stats[i].GetMemoryUse(); }
 
  os << "    Total memory consumed for geometry optimisation:   "
         << totalMemory/1024 << " kByte" << G4endl;
  os << "    Total CPU time elapsed for geometry optimisation: " 
         << std::setprecision(4) << totalCpuTime << " seconds"
         << std::setprecision(6) << G4endl;
 
  //
  // First list: sort by total CPU time
  //
  std::sort( stats.begin(), stats.end(),
    [](const G4SmartVoxelStat& a, const G4SmartVoxelStat& b)
  {
    return a.GetTotalTime() > b.GetTotalTime();
  } );
         
  const G4int maxPrint = 20;
  G4int nPrint = std::min ( nStat, maxPrint );

  if (nPrint != 0)
  {
    os << "\n    Voxelisation: top CPU users:" << G4endl;
    os << "    Percent   Total CPU    System CPU       Memory  Volume\n"
       << "    -------   ----------   ----------     --------  ----------"
       << G4endl;
  }

  for(i=0; i<nPrint; ++i)
  {
    G4double total = stats[i].GetTotalTime();
    G4double system = stats[i].GetSysTime();
    G4double perc = 0.0;

    if (system < 0) { system = 0.0; }
    if ((total < 0) || (totalCpuTime < perMillion))
      { total = 0; }
    else
      { perc = total*100/totalCpuTime; }

    os << std::setprecision(2) 
           << std::setiosflags(std::ios::fixed|std::ios::right)
           << std::setw(11) << perc
           << std::setw(13) << total
           << std::setw(13) << system
           << std::setw(13) << (stats[i].GetMemoryUse()+512)/1024
           << "k " << std::setiosflags(std::ios::left)
           << stats[i].GetVolume()->GetName()
           << std::resetiosflags(std::ios::floatfield|std::ios::adjustfield)
           << std::setprecision(6)
           << G4endl;
  }
 
  //
  // Second list: sort by memory use
  //
  std::sort( stats.begin(), stats.end(),
    [](const G4SmartVoxelStat& a, const G4SmartVoxelStat& b)
  {
    return a.GetMemoryUse() > b.GetMemoryUse();
  } );
 
  if (nPrint != 0)
  {
    os << "\n    Voxelisation: top memory users:" << G4endl;
    os << "    Percent     Memory      Heads    Nodes   Pointers    Total CPU    Volume\n"
       << "    -------   --------     ------   ------   --------   ----------    ----------"
       << G4endl;
  }

  for(i=0; i<nPrint; ++i)
  {
    G4long memory = stats[i].GetMemoryUse();
    G4double totTime = stats[i].GetTotalTime();
    if (totTime < 0) { totTime = 0.0; }

    os << std::setprecision(2) 
       << std::setiosflags(std::ios::fixed|std::ios::right)
       << std::setw(11) << G4double(memory*100)/G4double(totalMemory)
       << std::setw(11) << memory/1024 << "k "
       << std::setw( 9) << stats[i].GetNumberHeads()
       << std::setw( 9) << stats[i].GetNumberNodes()
       << std::setw(11) << stats[i].GetNumberPointers()
       << std::setw(13) << totTime << "    "
       << std::setiosflags(std::ios::left)
       << stats[i].GetVolume()->GetName()
       << std::resetiosflags(std::ios::floatfield|std::ios::adjustfield)
       << std::setprecision(6)
       << G4endl;
  }
  os << "--------------------------------------------------------------------------------"
     << G4endl << G4endl;
}

// ***************************************************************************
// Check whether parallel optimisation was requested -- static (class) data.
// ***************************************************************************
//
G4bool G4GeometryManager::IsParallelOptimisationConfigured()
{
  return fOptimiseInParallelConfigured;
}

// ***************************************************************************
// Report whether parallel optimisation is done -- static (class) data.
// ***************************************************************************
//
G4bool G4GeometryManager::IsParallelOptimisationFinished()
{
  return fParallelVoxelOptimisationFinished;
}
