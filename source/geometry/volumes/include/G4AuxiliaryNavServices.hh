// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AuxiliaryNavServices.hh,v 1.3 2000-04-25 16:15:02 gcosmo Exp $
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
// - Created: Paul Kent, Aug 96

#ifndef G4AuxiliaryNavServices_hh
#define G4AuxiliaryNavServices_hh

#include "geomdefs.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#include "G4AffineTransform.hh"


class G4AuxiliaryNavServices
{

 public:  // with description

   static G4bool CheckPointOnSurface( const G4VSolid* sampleSolid, 
				      const G4ThreeVector& localPoint, 
				      const G4ThreeVector* globalDirection, 
				      const G4AffineTransform& sampleTransform,
				      const G4bool  locatedOnEdge);
     //
     // Is the track (Point, direction) inside the solid sampleSolid ? 
     // Returns true if we are going to enter the volume, which is the case if:
     //   - The point is inside
     //   - The point is on the surface and the direction points inside
     //     or along it.
     // Else returns false

 private:
 
   G4bool testOne();

};

#include "G4AuxiliaryNavServices.icc"

#endif
