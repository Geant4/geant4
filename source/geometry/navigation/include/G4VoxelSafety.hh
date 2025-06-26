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
// G4VoxelSafety
//
// Class description:
//
// Utility for isotropic safety in volumes containing only G4PVPlacement
// daughter volumes for which voxels have been constructed.
// Implementation extracted and modified/adapted from G4VoxelNavigation class.

// Author: John Apostolakis (CERN), 30 April 2010
// --------------------------------------------------------------------
#ifndef G4VOXELSAFETY_HH
#define G4VOXELSAFETY_HH 1

#include "geomdefs.hh"
#include "G4NavigationHistory.hh"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"

#include "G4BlockingList.hh"

#include <vector>  // Required for voxel handling & voxel stack

class G4SmartVoxelNode;
class G4SmartVoxelHeader;

/**
 * @brief G4VoxelSafety is an utility class for the handling isotropic safety
 * in volumes containing only G4PVPlacement daughter volumes for which voxels
 * have been constructed.
 */

class G4VoxelSafety
{
  public:

    /**
     * Constructor and Destructor.
     */
    G4VoxelSafety();
    ~G4VoxelSafety();

    /**
     * Calculates the isotropic distance to the nearest boundary from the
     * specified point in the local coordinate system. 
     * The localpoint utilised must be within the current volume.
     *  @param[in] localPoint Local point.
     *  @param[in] currentPhysical Current physical volume.
     *  @param[in] maxLength Maximum length beyond which volumes are not checked.
     *  @returns Isotropic distance of given point to closest surface.
     */
    G4double ComputeSafety( const G4ThreeVector& localPoint,
                            const G4VPhysicalVolume& currentPhysical, 
                                  G4double maxLength = DBL_MAX );

    /**
     * Verbosity control.
     *  @note If level>0 && G4VERBOSE, printout can occur.
     */
    inline G4int GetVerboseLevel() const { return fVerbose; } 
    inline void  SetVerboseLevel(G4int level) { fVerbose = level; } 

  protected:

    /**
     * Cycles through levels of headers to process each node level.
     *  @param[in] pHead Voxel header.
     *  @param[in] localPoint Local point.
     *  @param[in] maxLength Maximum length beyond which volumes are not checked.
     *  @param[in] currentPhysical Current volume (used for debug printout).
     *  @param[in] distUpperDepth Upper square distance from voxel.
     *  @param[in] previousMinSafety Minimum distance beyond which not to look.
     *  @returns Isotropic distance of the point to closest volume in all nodes.
     */
    G4double SafetyForVoxelHeader( const G4SmartVoxelHeader* pHead,
                                   const G4ThreeVector& localPoint,
                                         G4double maxLength,
                                   const G4VPhysicalVolume& currentPhysical,
                                         G4double distUpperDepth = 0.0,
                                         G4double previousMinSafety = DBL_MAX );

    /**
     * Calculates the safety for volumes included in current Voxel Node.
     *  @param[in] curVoxelNode Voxel node.
     *  @param[in] localPoint Local point.
     *  @returns Isotropic distance of given point to closest volume in node.
     */
    G4double SafetyForVoxelNode( const G4SmartVoxelNode* curVoxelNode,
                                 const G4ThreeVector& localPoint ); 

  private:

    // ---- BEGIN State - values used during computation of Safety ------------

    /** Blocked volumes */
    G4BlockingList fBlockList;

    /** Cached pointer to mother logical volume */
    G4LogicalVolume* fpMotherLogical = nullptr;

    // ---- BEGIN Voxel Stack information -------------------------------------

    /** Voxel depth
     * @note fVoxelDepth==0+ => fVoxelAxisStack(0+) contains axes of voxel
     *       fVoxelDepth==-1 -> not in voxel.
     */
    G4int fVoxelDepth = -1;

    /** Voxel axes */
    std::vector<EAxis> fVoxelAxisStack;

    /** No slices per voxel at each level */
    std::vector<G4int> fVoxelNoSlicesStack;

    /** Width of voxels at each level */
    std::vector<G4double> fVoxelSliceWidthStack; 

    /** Node no point is inside at each level */
    std::vector<G4int> fVoxelNodeNoStack;    

    /** Voxel headers at each level */
    std::vector<const G4SmartVoxelHeader*> fVoxelHeaderStack;

    // ----- END Voxel Stack information --------------------------------------

    G4bool fCheck = false;
    G4int fVerbose = 0;
    G4double kCarTolerance;
};

#endif
