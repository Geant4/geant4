// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBSbox.cc,v 1.1 1999-01-07 16:09:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Box builder implementation (KidBox)
// OC 060996

#include "G4NURBSbox.hh"

	G4NURBSbox::G4NURBSbox(G4double DX, G4double DY, G4double DZ)
		:
		G4NURBS	(
			2, 2,	// linear along U and V
			4, 5	// line with two 90 degrees folds along U
				// rectangle along V (3 folds)
			)
				// let's it Generate regular knots vector
				// (note we are calling the second constructor)
		{
		t_indCtrlPt i = 0;


		CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
		CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
		CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 
		CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 

		CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
		CP(mpCtrlPts[i++],-DX, DY, DZ, 1 ); 
		CP(mpCtrlPts[i++],-DX, DY,-DZ, 1 ); 
		CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 

		CP(mpCtrlPts[i++], DX,-DY, DZ, 1 ); 
		CP(mpCtrlPts[i++],-DX,-DY, DZ, 1 ); 
		CP(mpCtrlPts[i++],-DX,-DY,-DZ, 1 ); 
		CP(mpCtrlPts[i++], DX,-DY,-DZ, 1 ); 

		CP(mpCtrlPts[i++], DX,-DY, DZ, 1 ); 
		CP(mpCtrlPts[i++], DX,-DY, DZ, 1 ); 
		CP(mpCtrlPts[i++], DX,-DY,-DZ, 1 ); 
		CP(mpCtrlPts[i++], DX,-DY,-DZ, 1 ); 

		CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
		CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
		CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 
		CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 

		


		
		}

const char*	G4NURBSbox::Whoami() const
		{
		return "Box (as a folded piece)";
		}





