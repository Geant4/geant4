// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBShexahedron.cc,v 1.1 1999-01-07 16:09:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// hexahedron implementation
//  OC 17 9 96

#include "G4NURBShexahedron.hh"

	G4NURBShexahedron::G4NURBShexahedron(const G4ThreeVector c [8])
		:
		G4NURBS	(
			2, 2,	// linear along U and V
			5, 4	// Square x polyline
			)
				// let's it Generate regular knots vector
		{
		// we need to define control points, indeed


		CP(mpCtrlPts[To1d(0,0)] ,  c[0].x(), c[0].y(), c[0].z(), 1 );
		CP(mpCtrlPts[To1d(1,0)] ,  c[1].x(), c[1].y(), c[1].z(), 1 );
		CP(mpCtrlPts[To1d(2,0)] ,  c[1].x(), c[1].y(), c[1].z(), 1 );
		CP(mpCtrlPts[To1d(3,0)] ,  c[0].x(), c[0].y(), c[0].z(), 1 );
		CP(mpCtrlPts[To1d(4,0)] ,  c[0].x(), c[0].y(), c[0].z(), 1 );


		CP(mpCtrlPts[To1d(0,1)] ,  c[0].x(), c[0].y(), c[0].z(), 1 );
		CP(mpCtrlPts[To1d(1,1)] ,  c[1].x(), c[1].y(), c[1].z(), 1 );
		CP(mpCtrlPts[To1d(2,1)] ,  c[2].x(), c[2].y(), c[2].z(), 1 );
		CP(mpCtrlPts[To1d(3,1)] ,  c[3].x(), c[3].y(), c[3].z(), 1 );
		CP(mpCtrlPts[To1d(4,1)] ,  c[0].x(), c[0].y(), c[0].z(), 1 );


		CP(mpCtrlPts[To1d(0,2)] ,  c[4].x(), c[4].y(), c[4].z(), 1 );
		CP(mpCtrlPts[To1d(1,2)] ,  c[5].x(), c[5].y(), c[5].z(), 1 );
		CP(mpCtrlPts[To1d(2,2)] ,  c[6].x(), c[6].y(), c[6].z(), 1 );
		CP(mpCtrlPts[To1d(3,2)] ,  c[7].x(), c[7].y(), c[7].z(), 1 );
		CP(mpCtrlPts[To1d(4,2)] ,  c[4].x(), c[4].y(), c[4].z(), 1 );


		CP(mpCtrlPts[To1d(0,3)] ,  c[4].x(), c[4].y(), c[4].z(), 1 );
		CP(mpCtrlPts[To1d(1,3)] ,  c[5].x(), c[5].y(), c[5].z(), 1 );
		CP(mpCtrlPts[To1d(2,3)] ,  c[5].x(), c[5].y(), c[5].z(), 1 );
		CP(mpCtrlPts[To1d(3,3)] ,  c[4].x(), c[4].y(), c[4].z(), 1 );
		CP(mpCtrlPts[To1d(4,3)] ,  c[4].x(), c[4].y(), c[4].z(), 1 );



		
		}

const char*	G4NURBShexahedron::Whoami() const
		{
		return "Hexahedron";
		}





