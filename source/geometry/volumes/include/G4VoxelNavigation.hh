// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VoxelNavigation.hh,v 1.3 1999-02-17 17:31:05 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4VoxelNavigation: Utility for navigation in volumes
// containing only G4PVPlacement daughter volumes for which voxels
// have been constructed.                         Paul Kent Aug 96

#ifndef G4VOXELNAVIGATION_HH
#define G4VOXELNAVIGATION_HH

#include "globals.hh"
#include "G4NavigationHistory.hh"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"

#include "G4BlockingList.hh"

// Required for voxel handling & voxel stack
#include "G4SmartVoxelProxy.hh"
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelHeader.hh"
#include <rw/tvvector.h>
#include <rw/tpvector.h>

// Voxel stack depth maximum [no resizing]
const G4int kNavigatorVoxelStackMax = 3;

class G4VoxelNavigation
{
public:
	G4VoxelNavigation();
        ~G4VoxelNavigation();
	G4SmartVoxelNode* VoxelLocate(G4SmartVoxelHeader *pHead,
                           const G4ThreeVector &localPoint);
	G4bool LevelLocate(G4NavigationHistory& history,
			   const G4VPhysicalVolume *blockedVol,
                           const G4int blockedNum,
                           const G4ThreeVector &globalPoint,
                           const G4ThreeVector* globalDirection,
			   const G4bool pLocatedOnEdge, 
			   G4ThreeVector &localPoint);

	G4double ComputeStep(const G4ThreeVector &globalPoint,
			     const G4ThreeVector &globalDirection,
			     const G4double currentProposedStepLength,
			     G4double &newSafety,
			     G4NavigationHistory &history,
			     G4bool &validExitNormal,
			     G4ThreeVector &exitNormal,
			     G4bool &exiting,
			     G4bool &entering,
                             G4VPhysicalVolume *(*pBlockedPhysical),
                             G4int &blockedReplicaNo);

        G4double ComputeSafety(const G4ThreeVector &globalpoint,
			       const G4NavigationHistory &history,
    		               const G4double pMaxLength=DBL_MAX );
private:
	G4double ComputeVoxelSafety(const G4ThreeVector &localPoint) const;
	G4bool LocateNextVoxel(const G4ThreeVector &localPoint,
			const G4ThreeVector& localDirection,
		        const G4double currentStep);
	G4bool VoxelSubLevelSetup(const G4ThreeVector& pLoc,
			G4SmartVoxelHeader *pHeader,
			G4int pDepth);

    G4BlockingList fBList;	// Blocked volumes
//
//  BEGIN Voxel Stack information
//

// Note: fVoxelDepth==0+ => fVoxelAxisStack(0+) contains axes of voxel
//       fVoxelDepth==-1 -> not in voxel
    G4int fVoxelDepth;

// Voxel axes
    RWTValVector<EAxis> fVoxelAxisStack;

// No slices per voxel at each level
    RWTValVector<G4int> fVoxelNoSlicesStack;

// Width of voxels at each level 
    RWTValVector<G4double> fVoxelSliceWidthStack; 

// Node no point is inside at each level 
    RWTValVector<G4int> fVoxelNodeNoStack;	  
				
// Voxel headers at each level
    RWTPtrVector<G4SmartVoxelHeader> fVoxelHeaderStack;

// Node containing last located point
    G4SmartVoxelNode* fVoxelNode;
//
//  END Voxel Stack information
//

};

#include "G4AuxiliaryNavServices.hh"   // Needed for inline implementation

#include "G4VoxelNavigation.icc"

#endif




