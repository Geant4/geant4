// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBScylinder.hh,v 2.1 1998/07/12 02:59:11 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Cylinder builder prototype
// OC 090796

#include "G4NURBS.hh"

#ifndef __C_G4NURBScylinder__

#define __C_G4NURBScylinder__ 1 

class	G4NURBScylinder : public G4NURBS
	{
	public:	G4NURBScylinder(G4double R, G4double DZ);
virtual	const char*	Whoami() const;
	};

#endif
// end of __C_G4NURBScylinder__




