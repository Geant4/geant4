// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NURBStubesector.cc,v 1.2 1999-05-12 16:11:02 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Tubesector builder implementation
// OC 290896

#include "G4NURBStubesector.hh"

// for sqrt
//#include <math.h>
// cf cylinder

// for ostrstream
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif


	G4NURBStubesector::G4NURBStubesector(G4double r, G4double R, G4double DZ, G4double PHI1, G4double PHI2)
		:
		G4NURBS	(
			2, 3,	// linear along U, quadratic along V
			5, DecideNbrCtrlPts(PHI1, PHI2),	
				// rectangle along U,  required stuff along V
				// we must use a static function which
				// take the two angles because the
				// mother constructor is initialised
				// before everything
			Regular,	// the knot vector along U will be generated
			RegularRep	// circular like knot vector also
			)
		{	

		// check angles
		G4double deltaPHI = PHI2-PHI1;
		while (deltaPHI <= 0) { PHI2 += 2*M_PI; deltaPHI += 2*M_PI; };

		G4int f = (int)floor(deltaPHI / (M_PI_2));	//number of pi/2 arcs
		
		const G4double mr = (r+R)/2;

		const G4double cp1 = cos(PHI1);
		const G4double sp1 = sin(PHI1);
		const G4double cp2 = cos(PHI2);
		const G4double sp2 = sin(PHI2);

		
		// define control points
		CP(mpCtrlPts[ 0] ,  cp1*mr, sp1*mr,  0, 1 );
		CP(mpCtrlPts[ 1] ,  cp1*mr, sp1*mr,  0, 1 );
		CP(mpCtrlPts[ 2] ,  cp1*mr, sp1*mr,  0, 1 );
		CP(mpCtrlPts[ 3] ,  cp1*mr, sp1*mr,  0, 1 );
		CP(mpCtrlPts[ 4] ,  cp1*mr, sp1*mr,  0, 1 );

		CP(mpCtrlPts[ 5] ,  cp1*mr, sp1*mr,  0, 1 );
		CP(mpCtrlPts[ 6] ,  cp1*mr, sp1*mr,  0, 1 );
		CP(mpCtrlPts[ 7] ,  cp1*mr, sp1*mr,  0, 1 );
		CP(mpCtrlPts[ 8] ,  cp1*mr, sp1*mr,  0, 1 );
		CP(mpCtrlPts[ 9] ,  cp1*mr, sp1*mr,  0, 1 );

		CP(mpCtrlPts[10] ,  cp1*r, sp1*r,  DZ, 1 );
		CP(mpCtrlPts[11] ,  cp1*R, sp1*R,  DZ, 1 );
		CP(mpCtrlPts[12] ,  cp1*R, sp1*R, -DZ, 1 );
		CP(mpCtrlPts[13] ,  cp1*r, sp1*r, -DZ, 1 );
		CP(mpCtrlPts[14] ,  cp1*r, sp1*r,  DZ, 1 );

		t_indCtrlPt	i = 15;
		G4double	srcAngle = PHI1;
		G4double	deltaAngleo2;

		G4double destAngle = M_PI_2 + PHI1;

		for(; f > 0; f--)
			{

			// the first arc CP is already Done

			deltaAngleo2 = (destAngle - srcAngle) / 2;
			const G4double csa = cos(srcAngle);
			const G4double ssa = sin(srcAngle);
			const G4double tdao2 = tan(deltaAngleo2); 

			// to calculate the intermediate CP :
			// rotate by srcAngle the (1, tdao2) point
			const t_Coord x = csa - ssa*tdao2;
			const t_Coord y = ssa + csa*tdao2;

			// weight of the CP
			const G4Float weight = (cos(deltaAngleo2));

			// initialization. postfix ++ because i initialized to 15
			CP(mpCtrlPts[i++], x*r, y*r,  DZ, 1, weight);
			CP(mpCtrlPts[i++], x*R, y*R,  DZ, 1, weight);
			CP(mpCtrlPts[i++], x*R, y*R, -DZ, 1, weight);
			CP(mpCtrlPts[i++], x*r, y*r, -DZ, 1, weight);
			CP(mpCtrlPts[i++], x*r, y*r,  DZ, 1, weight);

			// end CP (which is the first CP of the next arc)
			const G4double cda = cos(destAngle);
			const G4double sda = sin(destAngle);
			CP(mpCtrlPts[i++], cda*r, sda*r,  DZ, 1);
			CP(mpCtrlPts[i++], cda*R, sda*R,  DZ, 1);
			CP(mpCtrlPts[i++], cda*R, sda*R, -DZ, 1);
			CP(mpCtrlPts[i++], cda*r, sda*r, -DZ, 1);
			CP(mpCtrlPts[i++], cda*r, sda*r,  DZ, 1);

			// prepare next arc
			srcAngle = destAngle;
			destAngle += M_PI_2;
			};


			
		// f == 0, final Arc
		// could be handled in the loops

		destAngle = PHI2;
		deltaAngleo2 = (destAngle - srcAngle) / 2;
		const G4double csa = cos(srcAngle);
		const G4double ssa = sin(srcAngle);
		const G4double tdao2 = tan(deltaAngleo2); 

		// to calculate the intermediate CP :
		// rotate by srcAngle the (1, tdao2) point
		const t_Coord x = csa - ssa*tdao2;
		const t_Coord y = ssa + csa*tdao2;

		// weight of the CP
		const G4Float weight = (cos(deltaAngleo2));

		// initialization.
		CP(mpCtrlPts[i++], x*r, y*r,  DZ, 1, weight);
		CP(mpCtrlPts[i++], x*R, y*R,  DZ, 1, weight);
		CP(mpCtrlPts[i++], x*R, y*R, -DZ, 1, weight);
		CP(mpCtrlPts[i++], x*r, y*r, -DZ, 1, weight);
		CP(mpCtrlPts[i++], x*r, y*r,  DZ, 1, weight);

		// end CP
		const G4double cda = cos(destAngle);
		const G4double sda = sin(destAngle);
		CP(mpCtrlPts[i++], cda*r, sda*r,  DZ, 1);
		CP(mpCtrlPts[i++], cda*R, sda*R,  DZ, 1);
		CP(mpCtrlPts[i++], cda*R, sda*R, -DZ, 1);
		CP(mpCtrlPts[i++], cda*r, sda*r, -DZ, 1);
		CP(mpCtrlPts[i++], cda*r, sda*r,  DZ, 1);


/**/		if (i != (mtotnbrCtrlPts - 10) ) 
			{ G4cerr 
			<< "\nERROR: G4NURBStubesector::G4NURBStubesector: wrong index,"
			<< i << " instead of " << (mtotnbrCtrlPts - 10)
			<< "\n\tIt sounds very strange. The tubesector won't be correct. Have a nice debuging!" 
			<< endl;
			};

		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
		

		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
		CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);

		// possible to put a DZ DZ -DZ -DZ DZ column to scratch to a line instead of a point

		// creating the nurbs identity
		mpwhoami = new char [200];
		ostrstream	tmpstr(mpwhoami, 200);
		tmpstr << "Tubs" << " \tPHI1=" << PHI1 << " ; PHI2=" << PHI2 << '\0';
		// could be more sophisticated, reallocating
		// mpwhoami to the exact length
		
		}


const char*	G4NURBStubesector::Whoami() const
		{
		return mpwhoami;
		}

		G4NURBStubesector::~G4NURBStubesector()
		{
		if (mpwhoami) { delete [] mpwhoami; mpwhoami = NULL; };
		}


G4NURBStubesector::t_inddCtrlPt	G4NURBStubesector::DecideNbrCtrlPts(G4double PHI1, G4double PHI2)
		{
		// check angles
		G4double deltaPHI = PHI2-PHI1;
		while (deltaPHI <= 0) { PHI2 += 2*M_PI; deltaPHI += 2*M_PI; };
		G4double k = deltaPHI / (M_PI_2);

//		G4cerr << " k " << k << endl;
//		G4cerr << " fk " << floor(k) << endl;
//		G4cerr <<  " ifk " << ((int)(floor(k))) << endl;
//		G4cerr << " n " << (2*((int)(floor(k))) + 7) << endl;

		return ( 2*((int)(floor(k))) + 7 ); 		
		}
