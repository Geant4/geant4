// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBStube.hh,v 1.4.2.1 1999/12/07 20:48:49 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Tube builder prototype
// OC 090796

// Class Description:
// Tube builder prototype for NURBS.
// See documentation in graphics_reps/doc for details.
// Class Description - End:

#include "G4NURBS.hh"

#ifndef __C_G4NURBStube__

#define __C_G4NURBStube__ 1 

class	G4NURBStube : public G4NURBS
	{
	public:
	  G4NURBStube(G4double RMIN, G4double RMAX, G4double DZ);
	  virtual G4Visible&  operator = (const G4Visible& right);
	  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
	  const char*	Whoami() const;
	};

#endif
// end of __C_G4NURBStube__
