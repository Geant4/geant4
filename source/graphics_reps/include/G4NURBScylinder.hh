// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBScylinder.hh,v 1.3 1999-05-19 08:33:39 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
	public:
	  G4NURBScylinder(G4double R, G4double DZ);
	  virtual G4Visible&  operator = (const G4Visible& right);
	  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
	  const char*	Whoami() const;
	};

#endif
// end of __C_G4NURBScylinder__
