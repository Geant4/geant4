//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ParameterisedNavigation.hh,v 1.6 2002-07-23 08:50:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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
// ********************************************************************

#ifndef G4PARAMETERISEDNAVIGATION_HH
#define G4PARAMETERISEDNAVIGATION_HH

#include "globals.hh"

#include "g4std/vector"

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

    G4SmartVoxelNode* ParamVoxelLocate( G4SmartVoxelHeader* pHead,
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

    //  Voxel Stack information (for 1D optimisation only)
    //
    EAxis fVoxelAxis;
    G4int fVoxelNoSlices;
    G4double fVoxelSliceWidth; 
    G4int fVoxelNodeNo;  
    G4SmartVoxelHeader* fVoxelHeader;
};

#include "G4ParameterisedNavigation.icc"

#endif
