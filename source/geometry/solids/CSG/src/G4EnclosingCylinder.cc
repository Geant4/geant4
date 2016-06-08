//
// G4EnclosingCylinder.cc
//
// Implementation of a utility class for a quick check of geometry
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "G4EnclosingCylinder.hh"
#include "G4ReduciblePolygon.hh"

//
// Constructor
//
G4EnclosingCylinder::G4EnclosingCylinder( const G4ReduciblePolygon *rz,
					  const G4bool thePhiIsOpen, 
					  const G4double theStartPhi, const G4double theTotalPhi )
{
	//
	// Obtain largest r and smallest and largest z
	//
	radius = rz->Amax();
	zHi = rz->Bmax();
	zLo = rz->Bmin();
	
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
		
		concave = totalPhi > M_PI;
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
G4bool G4EnclosingCylinder::MustBeOutside( const G4ThreeVector &p ) const
{
	if (p.perp() > radius) return true;
	if (p.z() < zLo) return true;
	if (p.z() > zHi) return true;

	if (phiIsOpen) {
		if (concave) {
			if ( ((p.x()-dx1)*ry1 - (p.y()-dy1)*rx1) < 0) return false;
			if ( ((p.x()-dx2)*ry2 - (p.y()-dy2)*rx2) > 0) return false;
		}
		else {
			if ( ((p.x()-dx1)*ry1 - (p.y()-dy1)*rx1) > 0) return true;
			if ( ((p.x()-dx2)*ry2 - (p.y()-dy2)*rx2) < 0) return true;
		}
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
G4bool G4EnclosingCylinder::ShouldMiss( const G4ThreeVector &p, const G4ThreeVector &v ) const
{
	if (!MustBeOutside(p)) return false;
	
	G4double cross = p.x()*v.y() - p.y()*v.x();
	if (cross > radius) return true;
	
	if (p.perp() > radius) {
		G4double dot = p.x()*v.x() + p.y()*v.y();
		if (dot > 0) return true;
	}

	return false;
}		
