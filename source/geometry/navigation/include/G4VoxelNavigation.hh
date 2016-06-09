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
// $Id: G4VoxelNavigation.hh,v 1.5 2007/05/11 13:43:59 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// 
// class G4VoxelNavigation
//
// Class description:
//
// Utility for navigation in volumes containing only G4PVPlacement
// daughter volumes for which voxels have been constructed.

// History:
// - Created. Paul Kent, Aug 96
// --------------------------------------------------------------------
#ifndef G4VOXELNAVIGATION_HH
#define G4VOXELNAVIGATION_HH

#include "geomdefs.hh"
#include "G4NavigationHistory.hh"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"

#include "G4BlockingList.hh"

// Required for inline implementation
//
#include "G4AuxiliaryNavServices.hh"

// Required for voxel handling & voxel stack
//
#include <vector>
#include "G4SmartVoxelProxy.hh"
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelHeader.hh"

class G4VoxelNavigation
{
  public:  // with description

    G4VoxelNavigation();
    virtual ~G4VoxelNavigation();

    G4SmartVoxelNode* VoxelLocate( G4SmartVoxelHeader* pHead,
                             const G4ThreeVector& localPoint );

    virtual G4bool LevelLocate( G4NavigationHistory& history,
                          const G4VPhysicalVolume* blockedVol,
                          const G4int blockedNum,
                          const G4ThreeVector& globalPoint,
                          const G4ThreeVector* globalDirection,
                          const G4bool pLocatedOnEdge, 
                                G4ThreeVector& localPoint );

    virtual G4double ComputeStep( const G4ThreeVector& globalPoint,
                                  const G4ThreeVector& globalDirection,
                                  const G4double currentProposedStepLength,
                                        G4double& newSafety,
                                        G4NavigationHistory& history,
                                        G4bool& validExitNormal,
                                        G4ThreeVector& exitNormal,
                                        G4bool& exiting,
                                        G4bool& entering,
                                        G4VPhysicalVolume *(*pBlockedPhysical),
                                        G4int& blockedReplicaNo );

    virtual G4double ComputeSafety( const G4ThreeVector& globalpoint,
                                    const G4NavigationHistory& history,
                                    const G4double pMaxLength=DBL_MAX );

    inline G4int GetVerboseLevel() const;
    inline void  SetVerboseLevel(G4int level);
      // Get/Set Verbose(ness) level.
      // [if level>0 && G4VERBOSE, printout can occur]

    inline void  CheckMode(G4bool mode);
      // Run navigation in "check-mode", therefore using additional
      // verifications and more strict correctness conditions.
      // Is effective only with G4VERBOSE set.

  protected:

    G4double ComputeVoxelSafety( const G4ThreeVector& localPoint ) const;
    G4bool LocateNextVoxel( const G4ThreeVector& localPoint,
                            const G4ThreeVector& localDirection,
                            const G4double currentStep );

    G4BlockingList fBList;
      // Blocked volumes

    //
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

    std::vector<G4SmartVoxelHeader*> fVoxelHeaderStack;
      // Voxel headers at each level

    G4SmartVoxelNode* fVoxelNode;
      // Node containing last located point

    //
    //  END Voxel Stack information
    //

    G4bool fCheck;
    G4int  fVerbose;
    G4double kCarTolerance;
};

#include "G4VoxelNavigation.icc"

#endif
