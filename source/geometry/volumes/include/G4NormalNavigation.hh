// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NormalNavigation.hh,v 1.3 2000-04-25 16:15:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4NormalNavigation
//
// Class description:
//
// Utility for navigation in volumes containing only G4PVPlacement
// daughter volumes.

// History:
// - Created. Paul Kent, Aug 96

#ifndef G4NORMALNAVIGATION_HH
#define G4NORMALNAVIGATION_HH

#include "geomdefs.hh"
#include "G4NavigationHistory.hh"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"

class G4NormalNavigation
{
  public:  // with description

    G4bool LevelLocate(G4NavigationHistory &history,
		       const G4VPhysicalVolume *blockedVol,
                       const G4int blockedNum,
                       const G4ThreeVector &globalPoint,
                       const G4ThreeVector* globalDirection,
		       const G4bool pLocatedOnEdge, 
		       G4ThreeVector &localPoint);
      // Search positioned volumes in mother at current top level of history
      // for volume containing globalPoint. Do not test the blocked volume.
      // If a containing volume is found, `stack' the new volume and return
      // true, else return false (the point lying in the mother but not any
      // of the daughters). localPoint = global point in local system on entry,
      // point in new system on exit.

    G4double ComputeStep(const G4ThreeVector &localPoint,
			 const G4ThreeVector &localDirection,
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
};

#include "G4AuxiliaryNavServices.hh"    // Needed for inline methods

#include "G4NormalNavigation.icc"

#endif
