// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EnclosingCylinder.hh,v 1.1 2000-04-07 10:55:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4EnclosingCylinder
//
// Class description:
//
//   Definition of a utility class for quickly deciding if a point
//   is clearly outside a polyhedra or polycone or deciding if
//   a trajectory is clearly going to miss those shapes.

// --------------------------------------------------------------------

#ifndef G4EnclosingCylinder_hh
#define G4EnclosingCylinder_hh

#include "globals.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4ReduciblePolygon;

class G4EnclosingCylinder {
	public:
	G4EnclosingCylinder( const G4ReduciblePolygon *rz,
			     const G4bool phiIsOpen, 
			     const G4double startPhi, const G4double totalPhi );
	~G4EnclosingCylinder();
	
	G4bool MustBeOutside( const G4ThreeVector &p ) const;
	G4bool ShouldMiss( const G4ThreeVector &p, const G4ThreeVector &v ) const;
	
	protected:
	G4double radius; 	// radius of our cylinder
	G4double zLo, zHi;	// z extent
	
	G4bool 	 phiIsOpen;	// true if there is a phi segment
	G4double startPhi,	// for isPhiOpen==true, starting of phi segment
		 totalPhi;	// for isPhiOpen==true, size of phi segment

	G4double rx1, ry1,
		 dx1, dy1;
	G4double rx2, ry2,
		 dx2, dy2;
		 
	G4bool	 concave;	// True, if x/y cross section is concave
		
};

#endif
