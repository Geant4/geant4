//
// G4EnclosingCylinder.cc
//
// Implementation of a utility class for a quick check of geometry
//

#include "G4EnclosingCylinder.hh"

//
// Constructor
//
G4EnclosingCylinder::G4EnclosingCylinder( const G4double r[], const G4double z[], const  G4int n,
					  const G4bool thePhiIsOpen, 
					  const G4double theStartPhi, const G4double theTotalPhi )
{
	//
	// Obtain larges r and smallest and larges z
	//
	radius = r[0];
	zLo = zHi = z[0];
	const G4double *rr = r, *zz = z;
	while( ++zz, ++rr < r+n ) {
		if (*rr > radius) radius = *rr;
		if (*zz > zHi ) zHi = *zz;
		if (*zz < zLo ) zLo = *zz;
	}
	
	//
	// Save phi info
	//
	if ( phiIsOpen = thePhiIsOpen ) {
		startPhi = theStartPhi;
		totalPhi = theTotalPhi;
		
		rx1 = cos(startPhi);
		ry1 = sin(startPhi);
		dx1 = +ry1*10*kCarTolerance;
		dy1 = -rx1*10*kCarTolerance;
		
		rx2 = cos(startPhi+totalPhi);
		ry2 = sin(startPhi+totalPhi);
		dx2 = -ry2*10*kCarTolerance;
		dy2 = +rx2*10*kCarTolerance;
	}
	
	//
	// Add safety
	//
	radius += 10*kCarTolerance;
	zLo    -= 10*kCarTolerance;
	zHi    += 10*kCarTolerance;
}

//
// Destructor
//
G4EnclosingCylinder::~G4EnclosingCylinder() {;}


//
// Outside
//
// Decide very rapidly if the point is outside the cylinder
//
// If one is not certain, return false
//
G4bool G4EnclosingCylinder::Outside( const G4ThreeVector &p ) const
{
	if (p.perp() > radius) return true;
	if (p.z() < zLo) return true;
	if (p.z() > zHi) return true;

	if (phiIsOpen) {
		if ( ((p.x()-dx1)*ry1 - (p.y()-dy1)*rx1) > 0) return false;
		if ( ((p.x()-dx2)*ry2 - (p.y()-dy2)*rx2) < 0) return false;
		return true;
	}
	
	return false;
}
		
		
//
// Misses
//
// Decide very rapidly if the trajectory is going to miss the cylinder
//
// If one is not sure, return false
//
G4bool G4EnclosingCylinder::Misses( const G4ThreeVector &p, const G4ThreeVector &v ) const
{
	if (!Outside(p)) return false;
	
	G4double cross = p.x()*v.y() - p.y()*v.x();
	if (cross > radius) return true;
	
	if (p.perp() > radius) {
		G4double dot = p.x()*v.x() + p.y()*v.y();
		if (dot > 0) return true;
	}

	return false;
}		
