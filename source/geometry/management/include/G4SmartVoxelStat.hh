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
// $Id: G4SmartVoxelStat.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4SmartVoxelStat
//
// Class description:
//
// Stores information on the performance of the smart voxel algorithm
// for an individual logical volume.

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4SmartVoxelStat_hh
#define G4SmartVoxelStat_hh

#include "G4Types.hh"
#include <functional>

class G4LogicalVolume;
class G4SmartVoxelHeader;

class G4SmartVoxelStat
{
  public:  // with description
  
    G4SmartVoxelStat( const G4LogicalVolume *theVolume,
                      const G4SmartVoxelHeader *theVoxel,
                            G4double theSysTime,
                            G4double theUserTime );
      // Construct information on one volume's voxels

    const G4LogicalVolume *GetVolume() const;
      // Return a pointer to the logical volume
  
    const G4SmartVoxelHeader *GetVoxel() const;
      // Return a pointer to the voxel header
  
    G4double GetSysTime() const;
      // Get amount of system CPU time needed to build voxels
  
    G4double GetUserTime() const;
      // Get amount of user CPU time needed to build voxels
  
    G4double GetTotalTime() const;
      // Get total amount of CPU time needed to build voxels
  
    G4long GetNumberHeads() const;
      // Get number of voxel headers used in the volume
  
    G4long GetNumberNodes() const;
      // Get number of voxel slices used in the volume
  
    G4long GetNumberPointers() const;
      // Get number of voxel proxy pointers used in the volume
  
    G4long GetMemoryUse() const;
      // Get number of bytes needed to store voxel information


  protected:
  
    void CountHeadsAndNodes( const G4SmartVoxelHeader *head );
  
    const G4LogicalVolume *volume;
    const G4SmartVoxelHeader *voxel;
  
    G4double sysTime;
    G4double userTime;
  
    G4long heads;
    G4long nodes;
    G4long pointers;
  
  
  public:
  
    //
    // Functor objects for sorting
    //
    struct ByCpu
      : public std::binary_function< const G4SmartVoxelStat,
                                     const G4SmartVoxelStat, G4bool >
    {
      G4bool operator()( const G4SmartVoxelStat &a, const G4SmartVoxelStat &b )
      {
        return a.GetTotalTime() > b.GetTotalTime();
      }
    };
  
    struct ByMemory
      : public std::binary_function< const G4SmartVoxelStat,
                                     const G4SmartVoxelStat, G4bool >
    {
      G4bool operator()( const G4SmartVoxelStat &a, const G4SmartVoxelStat &b )
      {
        return a.GetMemoryUse() > b.GetMemoryUse();
      }
    };
};

#endif
