// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Parabola.cc,v 1.3 2000-08-28 08:57:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Parabola.cc
//
// ----------------------------------------------------------------------

#include "G4Parabola.hh"
#include "G4CurvePoint.hh"

G4Parabola::G4Parabola(){}
G4Parabola::~G4Parabola(){}

G4Curve* G4Parabola::Project(const G4Transform3D& tr)
{
  G4double axisZ= (tr*position.GetPZ()).unit().z();

  if (abs(axisZ)<kAngTolerance) 
    return 0;
  
  
  G4Vector3D newAxis(0, 0, axisZ>0? +1: -1);

  G4Vector3D xPrime= tr*position.GetPX();
  xPrime.setZ(0);
  G4Vector3D yPrime= tr*position.GetPY();
  yPrime.setZ(0);
  G4double u= -(xPrime*yPrime)/xPrime.mag2();

  G4Point3D newLocation= tr*position.GetLocation()+
    focalDist*(u*u*xPrime+2*u*yPrime);
  
  newLocation.setZ(0);
  G4Vector3D newRefDirection= xPrime;
  G4double newFocalDist= (focalDist*((2*u+1)*xPrime+2*yPrime)).mag()/sqrt(5.);
			  
  // create the new parabola
  G4Axis2Placement3D newPosition;
  newPosition.Init(newRefDirection, newAxis, newLocation);
  G4Parabola* r= new G4Parabola;
  r->Init(newPosition, newFocalDist);

  // introduce the shift in the parametrization
  // maybe the Sign must be changed?
  r->SetPShift(u);

  // set the bounds when necessary
  if (IsBounded()) 
    r->SetBounds(GetPStart(), GetPEnd());
  
  return r;
}


void G4Parabola::InitBounded()
{
  // the bbox must include the start and endpoints as well as the
  // extreme points if they lie on the curve
  bBox.Init(GetStart(), GetEnd());

  // the parameter values
  // belonging to the points with an extreme x, y and z coordinate
  for (G4int i=0; i<3; i++) 
  {
    G4double x_i= position.GetPX()(i);
    
    if (abs(x_i) <= kAngTolerance) 
    {
      G4double u= - position.GetPY()(i) / x_i;
      if (IsPOn(u)) 
	bBox.Extend(GetPoint(u));
    }
  }
}


G4bool G4Parabola::Tangent(G4CurvePoint& cp, G4Vector3D& v)
{
  // The tangent is computed from the 3D point representation
  // for all conics. An alternaive implementation (based on
  // the parametric point) might be worthwhile adding
  // for efficiency.
  
  const G4Axis2Placement3D& pos= *(GetPosition());
  G4Point3D p= pos.GetToPlacementCoordinates() * cp.GetPoint();

  v= p.y()*pos.GetPX() + (2*focalDist)*pos.GetPY();
  return true;
}
