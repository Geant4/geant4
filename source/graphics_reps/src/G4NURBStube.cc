// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBStube.cc,v 1.2 1999-05-19 08:33:48 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Tube builder implementation
// OC 090796

#include "G4NURBStube.hh"

// for sqrt
//#include <math.h>
// cf cylinder

	G4NURBStube::G4NURBStube(G4double r, G4double R, G4double DZ)
		:
		G4NURBS	(
			2, 3,	// linear along U, quadratic along V
			5, 9,	// rectangle along U, circle along V
			Regular,	// the knot vector along U will be generated
			RegularRep	// knot vector for the circle also
			)
		{
		// define control points

		const G4double sr2o2 = sqrt(2.)/2.;

		CP(mpCtrlPts[ 0] ,  r, 0,  DZ, 1 );
		CP(mpCtrlPts[ 1] ,  R, 0,  DZ, 1 );
		CP(mpCtrlPts[ 2] ,  R, 0, -DZ, 1 );
		CP(mpCtrlPts[ 3] ,  r, 0, -DZ, 1 );
		CP(mpCtrlPts[ 4] ,  r, 0,  DZ, 1 );

		CP(mpCtrlPts[ 5] ,  r, r,  DZ, 1 , sr2o2);
		CP(mpCtrlPts[ 6] ,  R, R,  DZ, 1 , sr2o2);
		CP(mpCtrlPts[ 7] ,  R, R, -DZ, 1 , sr2o2);
		CP(mpCtrlPts[ 8] ,  r, r, -DZ, 1 , sr2o2);
		CP(mpCtrlPts[ 9] ,  r, r,  DZ, 1 , sr2o2);

		CP(mpCtrlPts[10] ,  0, r,  DZ, 1 );
		CP(mpCtrlPts[11] ,  0, R,  DZ, 1 );
		CP(mpCtrlPts[12] ,  0, R, -DZ, 1 );
		CP(mpCtrlPts[13] ,  0, r, -DZ, 1 );
		CP(mpCtrlPts[14] ,  0, r,  DZ, 1 );

		CP(mpCtrlPts[15] , -r, r,  DZ, 1 , sr2o2);
		CP(mpCtrlPts[16] , -R, R,  DZ, 1 , sr2o2);
		CP(mpCtrlPts[17] , -R, R, -DZ, 1 , sr2o2);
		CP(mpCtrlPts[18] , -r, r, -DZ, 1 , sr2o2);
		CP(mpCtrlPts[19] , -r, r,  DZ, 1 , sr2o2);

		CP(mpCtrlPts[20] , -r, 0,  DZ, 1 );
		CP(mpCtrlPts[21] , -R, 0,  DZ, 1 );
		CP(mpCtrlPts[22] , -R, 0, -DZ, 1 );
		CP(mpCtrlPts[23] , -r, 0, -DZ, 1 );
		CP(mpCtrlPts[24] , -r, 0,  DZ, 1 );

		CP(mpCtrlPts[25] , -r,-r,  DZ, 1 , sr2o2);
		CP(mpCtrlPts[26] , -R,-R,  DZ, 1 , sr2o2);
		CP(mpCtrlPts[27] , -R,-R, -DZ, 1 , sr2o2);
		CP(mpCtrlPts[28] , -r,-r, -DZ, 1 , sr2o2);
		CP(mpCtrlPts[29] , -r,-R,  DZ, 1 , sr2o2);

		CP(mpCtrlPts[30] ,  0,-r,  DZ, 1 );
		CP(mpCtrlPts[31] ,  0,-R,  DZ, 1 );
		CP(mpCtrlPts[32] ,  0,-R, -DZ, 1 );
		CP(mpCtrlPts[33] ,  0,-r, -DZ, 1 );
		CP(mpCtrlPts[34] ,  0,-r,  DZ, 1 );

		CP(mpCtrlPts[35] ,  r,-r,  DZ, 1 , sr2o2);
		CP(mpCtrlPts[36] ,  R,-R,  DZ, 1 , sr2o2);
		CP(mpCtrlPts[37] ,  R,-R, -DZ, 1 , sr2o2);
		CP(mpCtrlPts[38] ,  r,-r, -DZ, 1 , sr2o2);
		CP(mpCtrlPts[39] ,  r,-r,  DZ, 1 , sr2o2);

		CP(mpCtrlPts[40] ,  r, 0,  DZ, 1 );
		CP(mpCtrlPts[41] ,  R, 0,  DZ, 1 );
		CP(mpCtrlPts[42] ,  R, 0, -DZ, 1 );
		CP(mpCtrlPts[43] ,  r, 0, -DZ, 1 );
		CP(mpCtrlPts[44] ,  r, 0,  DZ, 1 );

		
		}


G4Visible & G4NURBStube::operator = (const G4Visible &right) {
  return G4Visible::operator = (right);
}

G4VVisPrim & G4NURBStube::operator = (const G4VVisPrim &right) {
  return G4VVisPrim::operator = (right);
}

const char*	G4NURBStube::Whoami() const
		{
		return "Tube";
		}




