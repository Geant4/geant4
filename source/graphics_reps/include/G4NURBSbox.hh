// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBSbox.hh,v 1.5 1999-12-15 14:50:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Box builder prototype
// OC 060996

// Class Description:
// Box builder prototype for NURBS.
// See documentation in graphics_reps/doc for details.
// Class Description - End:

#include "G4NURBS.hh"

#ifndef __C_G4NURBSbox__

#define __C_G4NURBSbox__ 1 

class	G4NURBSbox : public G4NURBS
	{
        public:
	  G4NURBSbox(G4double DX, G4double DY, G4double DZ);
	  virtual G4Visible&  operator = (const G4Visible& right);
	  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
	  const char*	Whoami() const;
	};

#endif
// end of __C_G4NURBSbox__
