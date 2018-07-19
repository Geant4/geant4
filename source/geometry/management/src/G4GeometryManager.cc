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
// $Id: G4GeometryManager.cc 103235 2017-03-22 15:53:48Z gcosmo $
//
// class G4GeometryManager
//
// Implementation
//
// Author:
// 26.07.95 P.Kent Initial version, including optimisation Build
// --------------------------------------------------------------------

#include <iomanip>
#include "G4Timer.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"

#ifdef  G4GEOMETRY_VOXELDEBUG
#include "G4ios.hh"
#endif

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

// ***************************************************************************
// Static class data
// ***************************************************************************
//
G4ThreadLocal G4GeometryManager* G4GeometryManager::fgInstance = 0;
G4ThreadLocal G4bool G4GeometryManager::fIsClosed = false;

// ***************************************************************************
// Constructor. Set the geometry to be open
// ***************************************************************************
//
G4GeometryManager::G4GeometryManager() 
{
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4GeometryManager::~G4GeometryManager()
{
  fgInstance = 0;
  fIsClosed = false;
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
  if (!fIsClosed)
  {
    if (pVolume)
    {
      BuildOptimisations(pOptimise, pVolume);
    }
    else
    {
      BuildOptimisations(pOptimise, verbose);
    }
    fIsClosed=true;
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
  if (fIsClosed)
  {
    if (pVolume)
    {
      DeleteOptimisations(pVolume);
    }
    else
    {
      DeleteOptimisations();
    }
    fIsClosed=false;
  }
}

// ***************************************************************************
// Returns status of geometry
// ***************************************************************************
//
G4bool G4GeometryManager::IsGeometryClosed()
{
  return fIsClosed;
}

// ***************************************************************************
// Returns the instance of the singleton.
// Creates it in case it's called for the first time.
// ***************************************************************************
//
G4GeometryManager* G4GeometryManager::GetInstance()
{
  if (!fgInstance)
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
// Creates optimisation info. Builds all voxels if allOpts=true
// otherwise it builds voxels only for replicated volumes.
// ***************************************************************************
//
void G4GeometryManager::BuildOptimisations(G4bool allOpts, G4bool verbose)
{
   G4Timer timer;
   G4Timer allTimer;
   std::vector<G4SmartVoxelStat> stats;
   if (verbose)  { allTimer.Start(); }

   G4LogicalVolumeStore* Store = G4LogicalVolumeStore::GetInstance();
   G4LogicalVolume* volume;
   G4SmartVoxelHeader* head;
 
   for (size_t n=0; n<Store->size(); n++)
   {
     if (verbose) timer.Start();
     volume=(*Store)[n];
     // For safety, check if there are any existing voxels and
     // delete before replacement
     //
     head = volume->GetVoxelHeader();
     delete head;
     volume->SetVoxelHeader(0);
     if (    ( (volume->IsToOptimise())
            && (volume->GetNoDaughters()>=kMinVoxelVolumesLevel1&&allOpts) )
          || ( (volume->GetNoDaughters()==1)
            && (volume->GetDaughter(0)->IsReplicated()==true)
            && (volume->GetDaughter(0)->GetRegularStructureId()!=1) ) ) 
     {
#ifdef G4GEOMETRY_VOXELDEBUG
       G4cout << "**** G4GeometryManager::BuildOptimisations" << G4endl
              << "     Examining logical volume name = "
              << volume->GetName() << G4endl;
#endif
       head = new G4SmartVoxelHeader(volume);
       if (head)
       {
         volume->SetVoxelHeader(head);
       }
       else
       {
         std::ostringstream message;
         message << "VoxelHeader allocation error." << G4endl
                 << "Allocation of new VoxelHeader" << G4endl
                 << "        for volume " << volume->GetName() << " failed.";
         G4Exception("G4GeometryManager::BuildOptimisations()", "GeomMgt0003",
                     FatalException, message);
       }
       if (verbose)
       {
         timer.Stop();
         stats.push_back( G4SmartVoxelStat( volume, head,
                                            timer.GetSystemElapsed(),
                                            timer.GetUserElapsed() ) );
       }
     }
     else
     {
       // Don't create voxels for this node
#ifdef G4GEOMETRY_VOXELDEBUG
       G4cout << "**** G4GeometryManager::BuildOptimisations" << G4endl
              << "     Skipping logical volume name = " << volume->GetName()
              << G4endl;
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
// Creates optimisation info for the specified volumes subtree.
// ***************************************************************************
//
void G4GeometryManager::BuildOptimisations(G4bool allOpts,
                                           G4VPhysicalVolume* pVolume)
{
   if (!pVolume) { return; }

   // Retrieve the mother logical volume, if not NULL,
   // otherwise apply global optimisation for the world volume
   //
   G4LogicalVolume* tVolume = pVolume->GetMotherLogical();
   if (!tVolume) { return BuildOptimisations(allOpts, false); }

   G4SmartVoxelHeader* head = tVolume->GetVoxelHeader();
   delete head;
   tVolume->SetVoxelHeader(0);
   if (    ( (tVolume->IsToOptimise())
          && (tVolume->GetNoDaughters()>=kMinVoxelVolumesLevel1&&allOpts) )
        || ( (tVolume->GetNoDaughters()==1)
          && (tVolume->GetDaughter(0)->IsReplicated()==true) ) ) 
   {
     head = new G4SmartVoxelHeader(tVolume);
     if (head)
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
     G4cout << "**** G4GeometryManager::BuildOptimisations" << G4endl
            << "     Skipping logical volume name = " << tVolume->GetName()
            << G4endl;
#endif
   }

   // Scan recursively the associated logical volume tree
   //
  tVolume = pVolume->GetLogicalVolume();
  if (tVolume->GetNoDaughters())
  {
    BuildOptimisations(allOpts, tVolume->GetDaughter(0));
  }
}

// ***************************************************************************
// Removes all optimisation info.
// Loops over all logical volumes, deleting non-null voxels pointers,
// ***************************************************************************
//
void G4GeometryManager::DeleteOptimisations()
{
  G4LogicalVolume* tVolume = 0;
  G4LogicalVolumeStore* Store = G4LogicalVolumeStore::GetInstance();
  for (size_t n=0; n<Store->size(); n++)
  {
    tVolume=(*Store)[n];
    delete tVolume->GetVoxelHeader();
    tVolume->SetVoxelHeader(0);
  }
}

// ***************************************************************************
// Removes optimisation info for the specified subtree.
// Scans recursively all daughter volumes, deleting non-null voxels pointers.
// ***************************************************************************
//
void G4GeometryManager::DeleteOptimisations(G4VPhysicalVolume* pVolume)
{
  if (!pVolume) { return; }

  // Retrieve the mother logical volume, if not NULL,
  // otherwise global deletion to world volume.
  //
  G4LogicalVolume* tVolume = pVolume->GetMotherLogical();
  if (!tVolume) { return DeleteOptimisations(); }
  delete tVolume->GetVoxelHeader();
  tVolume->SetVoxelHeader(0);

  // Scan recursively the associated logical volume tree
  //
  tVolume = pVolume->GetLogicalVolume();
  if (tVolume->GetNoDaughters())
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
  if (G4SolidStore::GetInstance()->size())
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
                                     G4double totalCpuTime )
{
  G4cout << "G4GeometryManager::ReportVoxelStats -- Voxel Statistics"
         << G4endl << G4endl;
 
  //
  // Get total memory use
  //
  G4int i, nStat = stats.size();
  G4long totalMemory = 0;
 
  for( i=0;i<nStat;++i )  { totalMemory += stats[i].GetMemoryUse(); }
 
  G4cout << "    Total memory consumed for geometry optimisation:   "
         << totalMemory/1024 << " kByte" << G4endl;
  G4cout << "    Total CPU time elapsed for geometry optimisation: " 
         << std::setprecision(2) << totalCpuTime << " seconds"
         << std::setprecision(6) << G4endl;
 
  //
  // First list: sort by total CPU time
  //
  std::sort( stats.begin(), stats.end(), G4SmartVoxelStat::ByCpu() );
         
  G4int nPrint = nStat > 10 ? 10 : nStat;

  if (nPrint)
  {
    G4cout << "\n    Voxelisation: top CPU users:" << G4endl;
    G4cout << "    Percent   Total CPU    System CPU       Memory  Volume\n"
           << "    -------   ----------   ----------     --------  ----------"
           << G4endl;
    //         12345678901.234567890123.234567890123.234567890123k .
  }

  for(i=0;i<nPrint;++i)
  {
    G4double total = stats[i].GetTotalTime();
    G4double system = stats[i].GetSysTime();
    G4double perc = 0.0;

    if (system < 0) { system = 0.0; }
    if ((total < 0) || (totalCpuTime < perMillion))
      { total = 0; }
    else
      { perc = total*100/totalCpuTime; }

    G4cout << std::setprecision(2) 
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
  std::sort( stats.begin(), stats.end(), G4SmartVoxelStat::ByMemory() );
 
  if (nPrint)
  {
    G4cout << "\n    Voxelisation: top memory users:" << G4endl;
    G4cout << "    Percent     Memory      Heads    Nodes   Pointers    Total CPU    Volume\n"
           << "    -------   --------     ------   ------   --------   ----------    ----------"
           << G4endl;
    //         12345678901.2345678901k .23456789.23456789.2345678901.234567890123.   .
  }

  for(i=0;i<nPrint;++i)
  {
    G4long memory = stats[i].GetMemoryUse();
    G4double totTime = stats[i].GetTotalTime();
    if (totTime < 0) { totTime = 0.0; }

    G4cout << std::setprecision(2) 
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
}
