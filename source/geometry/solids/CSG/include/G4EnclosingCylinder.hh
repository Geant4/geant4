//
// G4EnclosingCylinder.hh
//
// Definition of a utility class for quickly deciding if a point
// is clearly outside a polyhedra or ploycone or deciding if
// a trajectory is clearly going to miss said shapes.
//
#ifndef G4EnclosingCylinder_hh
#define G4EnclosingCylinder_hh

#include "globals.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4EnclosingCylinder {
	public:
	G4EnclosingCylinder( const G4double r[], const G4double z[], const G4int n,
			     const G4bool phiIsOpen, 
			     const G4double startPhi, const G4double totalPhi );
	~G4EnclosingCylinder();
	
	G4bool Outside( const G4ThreeVector &p ) const;
	G4bool Misses( const G4ThreeVector &p, const G4ThreeVector &v ) const;
	
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
		
};

#endif
