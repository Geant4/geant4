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
// $Id: G4VoxelSafety.hh 92576 2015-09-04 10:26:26Z gcosmo $
//
// class G4VoxelSafety
//
// Class description:
//
// Utility for isotropic safety in volumes containing only G4PVPlacement
// daughter volumes for which voxels have been constructed.

// History:
// - Created. John Apostolakis,  30 April 2010
// --------------------------------------------------------------------
#ifndef G4VOXELSAFETY_HH
#define G4VOXELSAFETY_HH

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

class G4VoxelSafety
{
  public:  // with description

    G4VoxelSafety();
    ~G4VoxelSafety();

    G4SmartVoxelNode* VoxelLocate( G4SmartVoxelHeader* pHead,
                             const G4ThreeVector& localPoint );

    G4double ComputeSafety( const G4ThreeVector&     localPoint,
                            const G4VPhysicalVolume& currentPhysical, 
                                  G4double           maxLength=DBL_MAX );

    inline G4int GetVerboseLevel() const { return fVerbose; } 
    inline void  SetVerboseLevel(G4int level) { fVerbose= level; } 
      //
      // If level>0 && G4VERBOSE, printout can occur

  protected:

    G4double  SafetyForVoxelHeader( const G4SmartVoxelHeader* pHead,
                                    const G4ThreeVector& localPoint,
                                    G4double maxLength,
                                    const G4VPhysicalVolume& currentPhysical,
                                    G4double  distUpperDepth = 0.0,
                                    G4double  previousMinSafety= DBL_MAX
                                   );

    G4double  SafetyForVoxelNode(   const G4SmartVoxelNode *curVoxelNode,
                              const G4ThreeVector& localPoint ); 

    G4SmartVoxelNode* VoxelLocateLight( G4SmartVoxelHeader* pHead,
                                        const G4ThreeVector& localPoint ) const;
  private:

    // BEGIN State - values used during computation of Safety 
    // 
    G4BlockingList fBlockList;
      // Blocked volumes

    G4LogicalVolume* fpMotherLogical;

    //  BEGIN Voxel Stack information
    //
    G4int fVoxelDepth;
      // Note: fVoxelDepth==0+ => fVoxelAxisStack(0+) contains axes of voxel
      //       fVoxelDepth==-1 -> not in voxel

    std::vector<EAxis> fVoxelAxisStack;
      // Voxel axes

    std::vector<G4int> fVoxelNoSlicesStack;
      // No slices per voxel at each level

    std::vector<G4double> fVoxelSliceWidthStack; 
      // Width of voxels at each level 

    std::vector<G4int> fVoxelNodeNoStack;    
      // Node no point is inside at each level 

    std::vector<const G4SmartVoxelHeader*> fVoxelHeaderStack;
      // Voxel headers at each level

    G4SmartVoxelNode* fVoxelNode;
      // Node containing last located point

    //
    //  END Voxel Stack information

    G4bool fCheck;
    G4int  fVerbose;
    G4double kCarTolerance;
};

#endif
