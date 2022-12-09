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
// class G4ParameterisedNavigation
//
// Class description:
//
// Utility for navigation in volumes containing a single G4PVParameterised
// volume for which voxels for the replicated volumes have been constructed.
// [Voxels MUST be along one axis only: NOT refined]

// History:
// - Created. Paul Kent, Aug 96
// --------------------------------------------------------------------
#ifndef G4PARAMETERISEDNAVIGATION_HH
#define G4PARAMETERISEDNAVIGATION_HH

#include "G4Types.hh"

#include <vector>

#include "G4VoxelNavigation.hh"
#include "G4NavigationHistory.hh"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4BlockingList.hh"

class G4ParameterisedNavigation : public G4VoxelNavigation
{
  public:  // with description

    G4ParameterisedNavigation();
    ~G4ParameterisedNavigation();

    inline G4SmartVoxelNode* ParamVoxelLocate( G4SmartVoxelHeader* pHead,
                                         const G4ThreeVector& localPoint );

    G4bool LevelLocate( G4NavigationHistory& history,
                  const G4VPhysicalVolume* blockedVol,
                  const G4int blockedNum,
                  const G4ThreeVector& globalPoint,
                  const G4ThreeVector* globalDirection,
                  const G4bool pLocatedOnEdge, 
                        G4ThreeVector& localPoint );

    G4double ComputeStep( const G4ThreeVector& globalPoint,
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

    G4double ComputeSafety( const G4ThreeVector& localPoint,
                            const G4NavigationHistory& history,
                            const G4double pProposedMaxLength=DBL_MAX );

  private:

    G4double ComputeVoxelSafety( const G4ThreeVector& localPoint,
                                 const EAxis pAxis ) const;
    G4bool LocateNextVoxel( const G4ThreeVector& localPoint,
                            const G4ThreeVector& localDirection,
                            const G4double currentStep,
                            const EAxis pAxis );
  private:

    // Necessary to resolve cases with nested parameterisations

    inline G4VSolid* IdentifyAndPlaceSolid( G4int num,
                                     G4VPhysicalVolume* apparentPhys, 
                                     G4VPVParameterisation* curParam );
       // Call virtual 'Compute' methods, and copy information if nested.
       // 'ApparentPhys' is potentially a PhysV or PhysT.

    G4VPhysicalVolume* CreateVolumeWithParent(G4VPhysicalVolume* curPhysical,
                                              const G4NavigationHistory& hist );
       // Create necessary parent touchable and physical with parent.

  private:

    //  Voxel Stack information (for 1D optimisation only)
    //
    EAxis fVoxelAxis = kUndefined;
    G4int fVoxelNoSlices = 0;
    G4double fVoxelSliceWidth = 0.0; 
    std::size_t fVoxelNodeNo = 0;  
    G4SmartVoxelHeader* fVoxelHeader = nullptr;
};

#include "G4ParameterisedNavigation.icc"

#endif
