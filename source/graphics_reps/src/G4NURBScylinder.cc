// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBScylinder.cc,v 1.3 1999-12-15 14:50:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Cylinder builder implementation
// OC 090796

#include "G4NURBScylinder.hh"

// for sqrt
//#include <math.h>
// math.h included in template.hh included in globals.hh included in
// G4NURBS.hh included in G4NURBScylinder.hh

// the cylinder constructor use the first G4NURBS constructor
// look in G4NURBStube if you want to see how to use the second one.

	G4NURBScylinder::G4NURBScylinder(G4double R, G4double DZ)
		:
		G4NURBS	(
			2, 3,	// linear along U, quadratic along V
			4, 9,	// half rectangle along U, circle along V
			(new t_CtrlPt [ 4 * 9 ]),  // the array for CtrlPts
			NULL,	// the knot vector along U will be generated
			(new t_Knot [ 3 + 9 ])   // knot vector for the circle
			)
		{
		// define the V knot vector
		m[V].pKnots[ 0] = 0;
		m[V].pKnots[ 1] = 0;
		m[V].pKnots[ 2] = 0;
		m[V].pKnots[ 3] = 0.25;
		m[V].pKnots[ 4] = 0.25;
		m[V].pKnots[ 5] = 0.5;
		m[V].pKnots[ 6] = 0.5;
		m[V].pKnots[ 7] = 0.75;
		m[V].pKnots[ 8] = 0.75;
		m[V].pKnots[ 9] = 1;
		m[V].pKnots[10] = 1;
		m[V].pKnots[11] = 1;

		// define control points

		const G4double sr2o2 = sqrt(2.)/2. ;

		CP(mpCtrlPts[ 0] ,  0,  0,  DZ, 1 );
		CP(mpCtrlPts[ 1] ,  R, 0,  DZ, 1 );
		CP(mpCtrlPts[ 2] ,  R, 0, -DZ, 1 );
		CP(mpCtrlPts[ 3] ,  0,  0, -DZ, 1 );

		CP(mpCtrlPts[ 4] ,  0,  0,  DZ, 1 );
		CP(mpCtrlPts[ 5] ,  R, R, DZ, 1 , sr2o2);
		CP(mpCtrlPts[ 6] ,  R, R,-DZ, 1 , sr2o2);
		CP(mpCtrlPts[ 7] ,  0,  0, -DZ, 1 );

		CP(mpCtrlPts[ 8] ,  0,  0,  DZ, 1 );
		CP(mpCtrlPts[ 9] ,  0, R,  DZ, 1 );
		CP(mpCtrlPts[10] ,  0, R, -DZ, 1 );
		CP(mpCtrlPts[11] ,  0,  0, -DZ, 1 );

		CP(mpCtrlPts[12] ,  0,  0,  DZ, 1 );
		CP(mpCtrlPts[13] , -R, R, DZ, 1 , sr2o2);
		CP(mpCtrlPts[14] , -R, R,-DZ, 1 , sr2o2);
		CP(mpCtrlPts[15] ,  0,  0, -DZ, 1 );

		CP(mpCtrlPts[16] ,  0,  0,  DZ, 1 );
		CP(mpCtrlPts[17] , -R, 0,  DZ, 1 );
		CP(mpCtrlPts[18] , -R, 0, -DZ, 1 );
		CP(mpCtrlPts[19] ,  0,  0, -DZ, 1 );

		CP(mpCtrlPts[20] ,  0,  0,  DZ, 1 );
		CP(mpCtrlPts[21] , -R,-R, DZ, 1 , sr2o2);
		CP(mpCtrlPts[22] , -R,-R,-DZ, 1 , sr2o2);
		CP(mpCtrlPts[23] ,  0,  0, -DZ, 1 );

		CP(mpCtrlPts[24] ,  0,  0,  DZ, 1 );
		CP(mpCtrlPts[25] ,  0, -R, DZ, 1 );
		CP(mpCtrlPts[26] ,  0, -R,-DZ, 1 );
		CP(mpCtrlPts[27] ,  0,  0, -DZ, 1 );

		CP(mpCtrlPts[28] ,  0,  0,  DZ, 1 );
		CP(mpCtrlPts[29] ,  R,-R, DZ, 1 , sr2o2);
		CP(mpCtrlPts[30] ,  R,-R,-DZ, 1 , sr2o2);
		CP(mpCtrlPts[31] ,  0,  0, -DZ, 1 );

		CP(mpCtrlPts[32] ,  0,  0,  DZ, 1 );
		CP(mpCtrlPts[33] ,  R, 0,  DZ, 1 );
		CP(mpCtrlPts[34] ,  R, 0, -DZ, 1 );
		CP(mpCtrlPts[35] ,  0,  0, -DZ, 1 );
		
		}

G4Visible & G4NURBScylinder::operator = (const G4Visible &right) {
  return G4Visible::operator = (right);
}

G4VVisPrim & G4NURBScylinder::operator = (const G4VVisPrim &right) {
  return G4VVisPrim::operator = (right);
}

const char*	G4NURBScylinder::Whoami() const
		{
		return "Cylinder";
		}






