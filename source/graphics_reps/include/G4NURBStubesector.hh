// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBStubesector.hh,v 1.2 1999-05-19 08:33:40 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Tubesector builder prototype
// OC 290896

#include "G4NURBS.hh"

#ifndef __C_G4NURBStubesector__

#define __C_G4NURBStubesector__ 1 

class	G4NURBStubesector : public G4NURBS
	{
	public:	
		// angle in radians

		// If PHI2 smaller (or equal) than PHI1 , it is incremented
		// by 2pi as necessary to become strictly greater.

		// Except that, you can use any value for the arguments,
		// it's the renderer or you that will have troubles. 

		G4NURBStubesector(G4double RMIN, G4double RMAX,
				  G4double DZ, G4double PHI1, G4double PHI2);
	  virtual 	~G4NURBStubesector();
	  virtual G4Visible&  operator = (const G4Visible& right);
	  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
	  virtual	const char*	Whoami() const;

	protected:
	char *	mpwhoami;

	private:
static	t_inddCtrlPt	DecideNbrCtrlPts(G4double PHI1, G4double PHI2);
	};




#endif
// end of __C_G4NURBStubesector__




