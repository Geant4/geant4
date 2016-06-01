// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBStube.hh,v 2.1 1998/07/12 02:59:13 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Tube builder prototype
// OC 090796

#include "G4NURBS.hh"

#ifndef __C_G4NURBStube__

#define __C_G4NURBStube__ 1 

class	G4NURBStube : public G4NURBS
	{
	public:	G4NURBStube(G4double RMIN, G4double RMAX, G4double DZ);
virtual	const char*	Whoami() const;
	};

#endif
// end of __C_G4NURBStube__




