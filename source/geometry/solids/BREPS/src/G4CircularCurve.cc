// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CircularCurve.cc,v 1.2 1999-12-15 14:50:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4CircularCurve.hh"
#include "G4Ellipse.hh"

// G4CircularCurve
G4CircularCurve::G4CircularCurve() {}
G4CircularCurve::~G4CircularCurve() {}

//////////////////////////////////////////////////////////////////////////////

void G4CircularCurve::InitBounded() 
{  
  // the bbox must include the start and endpoints as well as the
  // extreme points if they lie on the curve
  bBox.Init(GetStart(), GetEnd());
  
  // the parameter values
  // belonging to the points with an extreme x, y and z coordinate
  for (G4int i=0; i<3; i++) 
  {
    G4double u = atan2(position.GetPY()(i), position.GetPX()(i));
    
    if (IsPOn(u)) 
      bBox.Extend(GetPoint(u));
    
    if (IsPOn(u+pi)) 
      bBox.Extend(GetPoint(u+pi));
  }
}

//////////////////////////////////////////////////////////////////////////////

G4Curve* G4CircularCurve::Project(const G4Transform3D& tr)
{
  G4Ellipse e;
  e.Init(position, radius, radius);
  e.SetBounds(GetPStart(), GetPEnd());
 
  return e.Project(tr);
}

//////////////////////////////////////////////////////////////////////////////

G4bool G4CircularCurve::Tangent(G4CurvePoint& cp, G4Vector3D& v)
{
  // The tangent is computed from the 3D point representation
  // for all conics. An alternaive implementation (based on
  // the parametric point) might be worthwhile adding
  // for efficiency.
  
  const G4Axis2Placement3D& pos= *(GetPosition());
  G4Point3D p= pos.GetToPlacementCoordinates() * cp.GetPoint();
  
  v= -p.y()*pos.GetPX() + p.x()*pos.GetPY();
  return true;
}

