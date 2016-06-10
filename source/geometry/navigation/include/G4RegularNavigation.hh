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
// $Id: G4RegularNavigation.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// class G4RegularNavigation
//
// Class description:
//
// Utility for fast navigation in volumes containing a regular
// parameterisation. If two contiguous voxels have the same material,
// navigation does not stop at the surface

// History:
// - Created.   P. Arce, May 2007
// --------------------------------------------------------------------
#ifndef G4RegularNavigation_HH
#define G4RegularNavigation_HH

#include <vector>

#include "G4Types.hh"
#include "G4ThreeVector.hh"

class G4NormalNavigation;
class G4VPhysicalVolume;
class G4Navigator;
class G4NavigationHistory;

class G4RegularNavigation
{
  public:  // with description
  
    G4RegularNavigation();
   ~G4RegularNavigation();
  
    G4bool LevelLocate(      G4NavigationHistory& history,
                       const G4VPhysicalVolume* blockedVol,
                       const G4int blockedNum,
                       const G4ThreeVector& globalPoint,
                       const G4ThreeVector* globalDirection,
                       const G4bool pLocatedOnEdge, 
                             G4ThreeVector& localPoint );
      // Locate point using its position with respect to regular
      // parameterisation container volume.

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
      // Method never called because to be called the daughter has to be a
      // 'regular' volume. This would only happen if the track is in the
      // mother of voxels volume. But the voxels fill completely their mother,
      // so when a track enters the mother it automatically enters a voxel.
  
    G4double ComputeStepSkippingEqualMaterials( 
                                G4ThreeVector& localPoint,
                          const G4ThreeVector& globalDirection,
                          const G4double currentProposedStepLength,
                                G4double& newSafety,
                                G4NavigationHistory& history,
                                G4bool& validExitNormal,
                                G4ThreeVector& exitNormal,
                                G4bool& exiting,
                                G4bool& entering,
                                G4VPhysicalVolume *(*pBlockedPhysical),
                                G4int& blockedReplicaNo,
                                G4VPhysicalVolume* pCurrentPhysical);
      // Compute the step skipping surfaces when they separate voxels with
      // equal materials. Loop to voxels until a different material is found:
      // invokes G4NormalNavigation::ComputeStep() in each voxel and move the
      // point to the next voxel.

    G4double ComputeSafety( const G4ThreeVector& localPoint,
                            const G4NavigationHistory& history,
                            const G4double pProposedMaxLength=DBL_MAX );
      // Method never called because to be called the daughter has to be a
      // 'regular' volume. This would only happen if the track is in the
      // mother of voxels volume. But the voxels fill completely their mother,
      // so when a track enters the mother it automatically enters a voxel.

  public:  // without description

    // Set and Get methods

    void SetVerboseLevel(G4int level) { fverbose = level; }
    void CheckMode(G4bool mode) { fcheck = mode; }
    void SetNormalNavigation( G4NormalNavigation* fnormnav )
      { fnormalNav = fnormnav; }

  private:

    G4int fverbose;
    G4bool fcheck;

    G4NormalNavigation* fnormalNav;
    G4double kCarTolerance;  
};

#endif
