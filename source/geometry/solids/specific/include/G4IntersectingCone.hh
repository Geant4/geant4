// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IntersectingCone.hh,v 1.1 2000-04-07 10:56:27 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4IntersectingCone
//
// Class description:
//
//   Utility class which calculates the intersection
//   of an arbitrary line with a fixed cone
//
// --------------------------------------------------------------------

#ifndef G4IntersectingCone_hh
#define G4IntersectingCone_hh

#include "globals.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4IntersectingCone {
	public:
	G4IntersectingCone( const G4double r[2], const G4double z[2] );
	virtual ~G4IntersectingCone();
	
	G4int LineHitsCone( const G4ThreeVector &p, const G4ThreeVector &v,
                            G4double *s1, G4double *s2 );
	
	G4bool HitOn( const G4double r, const G4double z );
	
	inline G4double RLo() const { return rLo; }
	inline G4double RHi() const { return rHi; }
	inline G4double ZLo() const { return zLo; }
	inline G4double ZHi() const { return zHi; }
	
	
	protected:
	G4double zLo, zHi,	// Z bounds of side
		 rLo, rHi;	// R bounds of side

	G4bool	 type1;		// True if cone is type 1 (abs(z1-z2)>abs(r1-r2))
	G4double A, B;		// Cone radius parameter:
				//	type 1: r = A + B*z
				//	type 2: z = A + B*r

	G4int LineHitsCone1( const G4ThreeVector &p, const G4ThreeVector &v,
                             G4double *s1, G4double *s2 );
	G4int LineHitsCone2( const G4ThreeVector &p, const G4ThreeVector &v,
                             G4double *s1, G4double *s2 );
};

#endif
