// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Conic.cc,v 1.3 2000-08-28 08:57:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4conic.cc
//
// ----------------------------------------------------------------------

#include "G4Conic.hh"

G4Conic::G4Conic () : pShift(0)
{
}

G4Conic::~G4Conic()
{
}

/*
void
G4ConicalCurve::ProjectCurve(const G4Plane& Pl1, const G4Plane& Pl2)
{  
  // Curve start
  Project(ProjStart, Start, Pl1, Pl2);
  // Curve end
  Project(ProjEnd, End, Pl1, Pl2);
  // Placement
  Position.ProjectPlacement(Pl1,Pl2);
}

G4int
G4ConicalCurve::HitPartOfCurve(G4double Angle, G4double Solution,
                               const G4Point2d& ProjHit)
{
    	// Check if Solution1 is part of the curve i.e. in the "pie"
	G4double TmpSol1 = Solution - ProjHit.X();
	G4Point2d ArcHit1(TmpSol1, ProjHit.Y());
	G4double Cross1 = CrossProduct( ProjStart, ArcHit1);
	G4double Cross2 = CrossProduct( ArcHit1  , ProjEnd);
	if( (Angle<=0 && Cross1<=0 && Cross2 <=0) ||
	    (Angle> 0 && !(Cross1>=0 && Cross2 >=0)) )
	    // Solution1 is on the curve
	  return 1;
	return 0;
}
*/
