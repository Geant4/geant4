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
// G4SmartVoxelStat
//
// Class description:
//
// Stores information on the performance of the smart voxel algorithm
// for an individual logical volume.

// Author: D.C.Williams (UCSC), 1998
// --------------------------------------------------------------------
#ifndef G4SmartVoxelStat_hh
#define G4SmartVoxelStat_hh

#include "G4Types.hh"

class G4LogicalVolume;
class G4SmartVoxelHeader;

/**
 * @brief G4SmartVoxelStat stores the information on the performance of the
 * smart voxel optimisation algorithm for an individual logical volume.
 */

class G4SmartVoxelStat
{
  public:
  
    /**
     * Constructs the information on one volume's voxels.
     *  @param[in] theVolume Pointer to the logical volume concerned.
     *  @param[in] theVoxel Pointer to the associated voxel header.
     *  @param[in] theSysTime System time.
     *  @param[in] theUserTime User time.
     */
    G4SmartVoxelStat( const G4LogicalVolume* theVolume,
                      const G4SmartVoxelHeader* theVoxel,
                            G4double theSysTime,
                            G4double theUserTime );

    /**
     * Returns a pointer to the logical volume.
     */
    const G4LogicalVolume* GetVolume() const;
  
    /**
     * Returns a pointer to the voxel header.
     */
    const G4SmartVoxelHeader* GetVoxel() const;
  
    /**
     * Gets the amount of system CPU time needed to build voxels.
     */
    G4double GetSysTime() const;
  
    /**
     * Gets the amount of user CPU time needed to build voxels.
     */
    G4double GetUserTime() const;
  
    /**
     * Gets the total amount of CPU time needed to build voxels.
     */
    G4double GetTotalTime() const;
  
    /**
     * Gets the number of voxel headers used in the volume.
     */
    G4long GetNumberHeads() const;
  
    /**
     * Gets the number of voxel slices used in the volume.
     */
    G4long GetNumberNodes() const;

    /**
     * Gets the number of voxel proxy pointers used in the volume.
     */
    G4long GetNumberPointers() const;
  
    /**
     * Gets the number of bytes needed to store voxel information.
     */
    G4long GetMemoryUse() const;


  private:
  
    /**
     * Counts headers and nodes from provided 'head' header.
     */
    void CountHeadsAndNodes( const G4SmartVoxelHeader* head );
  
  private:

    const G4LogicalVolume* volume;
    const G4SmartVoxelHeader* voxel;
  
    G4double sysTime;
    G4double userTime;
  
    G4long heads = 1;
    G4long nodes = 0;
    G4long pointers = 0;
};

#endif
