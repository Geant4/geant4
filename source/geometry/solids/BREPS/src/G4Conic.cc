#include "G4Conic.hh"
// G4Conic

G4Conic::G4Conic (): pShift(0) {}
G4Conic::~G4Conic() {}



/*
void G4ConicalCurve::ProjectCurve(const G4Plane& Pl1, const G4Plane& Pl2)
{  
  // Curve start
  Project(ProjStart, Start, Pl1, Pl2);
  // Curve end
  Project(ProjEnd, End, Pl1, Pl2);
  // Placement
  Position.ProjectPlacement(Pl1,Pl2);
}

int G4ConicalCurve::HitPartOfCurve(G4double Angle, G4double Solution, const G4Point2d& ProjHit)
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
