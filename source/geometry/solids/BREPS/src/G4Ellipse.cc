// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Ellipse.cc,v 1.6 2000-11-10 17:41:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Ellipse.cc
//
// ----------------------------------------------------------------------

#include "G4Ellipse.hh"

G4Ellipse::G4Ellipse()
{
}

G4Ellipse::~G4Ellipse()
{
}

G4Ellipse::G4Ellipse(const G4Ellipse& right)
  : semiAxis1(right.semiAxis1), semiAxis2(right.semiAxis2),
    ratioAxis2Axis1(right.ratioAxis2Axis1), toUnitCircle(right.toUnitCircle),
    forTangent(right.forTangent)
{
  pShift    = right.pShift;
  position  = right.position;
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;
}

G4Ellipse& G4Ellipse::operator=(const G4Ellipse& right)
{
  if (&right == this) return *this;

  semiAxis1 = right.semiAxis1;
  semiAxis2 = right.semiAxis2;
  ratioAxis2Axis1 = right.ratioAxis2Axis1;
  toUnitCircle    = right.toUnitCircle;
  forTangent = right.forTangent;
  pShift   = right.pShift;
  position = right.position;
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;

  return *this;
}

G4Curve* G4Ellipse::Project(const G4Transform3D& tr)
{
  G4Point3D newLocation = tr*position.GetLocation();
  newLocation.setZ(0);
  G4double axisZ        = ( tr*position.GetPZ() ).unit().z();

  if (abs(axisZ)<kAngTolerance) 
    return 0;
  
  G4Vector3D newAxis(0, 0, axisZ>0? +1: -1);

  // get the parameter of an endpoint of an axis
  // (this is a point the distance of which from the center is extreme)
  G4Vector3D xPrime= tr*position.GetPX();
  xPrime.setZ(0);
  G4Vector3D yPrime= tr*position.GetPY();
  yPrime.setZ(0);

  G4Vector3D a = G4Vector3D( semiAxis1*xPrime );
  G4Vector3D b = G4Vector3D( semiAxis2*yPrime );
  
  G4double u;
  G4double abmag = a.mag2()-b.mag2();
  G4double prod = 2*a*b;

  if ((abmag > FLT_MAX) && (prod < -FLT_MAX))
    u = -pi/8;
  else if ((abmag < -FLT_MAX) && (prod > FLT_MAX))
    u = 3*pi/8;
  else if ((abmag < -FLT_MAX) && (prod < -FLT_MAX))
    u = -3*pi/8;
  else if ((abs(abmag) < perMillion) && (abs(prod) < perMillion))
    u = 0.;
  else
    u = atan2(prod,abmag) / 2;

  // get the coordinate axis directions and the semiaxis lengths
  G4Vector3D sAxis1          = G4Vector3D( a*cos(u)+b*sin(u) );
  G4Vector3D sAxis2          = G4Vector3D( a*cos(u+pi/2)+b*sin(u+pi/2) );
  G4double newSemiAxis1      = sAxis1.mag();
  G4double newSemiAxis2      = sAxis2.mag();
  G4Vector3D newRefDirection = sAxis1;

  // create the new ellipse
  G4Axis2Placement3D newPosition;
  newPosition.Init(newRefDirection, newAxis, newLocation);
  G4Ellipse* r= new G4Ellipse;
  r->Init(newPosition, newSemiAxis1, newSemiAxis2);
  
  // introduce the shift in the parametrization
  // maybe the Sign must be changed?
  r->SetPShift(u);
  
  // set the bounds when necessary
  if (IsBounded()) 
    r->SetBounds(GetPStart(), GetPEnd());
  
  // L. Broglia
  // copy sense of the curve
  r->SetSameSense(GetSameSense());

  return r;
}

void G4Ellipse::InitBounded()
{
  // original implementation
  // const G4Point3D& center = position.GetLocation();
  // G4double maxEntent      = G4std::max(semiAxis1, semiAxis2);
  // G4Vector3D halfExtent(maxEntent, maxEntent, maxEntent);
  // bBox.Init(center+halfExtent, center-halfExtent);

  // the bbox must include the start and endpoints as well as the
  // extreme points if they lie on the curve
  bBox.Init(GetStart(), GetEnd());

  // the parameter values
  // belonging to the points with an extreme x, y and z coordinate
  for (G4int i=0; i<3; i++) 
  {
    G4double u= atan2(position.GetPY()(i)*semiAxis2,
		      position.GetPX()(i)*semiAxis1);
    if (IsPOn(u)) 
      bBox.Extend(GetPoint(u));
    
    if (IsPOn(u+pi)) 
      bBox.Extend(GetPoint(u+pi));
  }
}


G4bool G4Ellipse::Tangent(G4CurvePoint& cp, G4Vector3D& v)
{
  // The tangent is computed from the 3D point representation
  // for all conics. An alternaive implementation (based on
  // the parametric point) might be worthwhile adding
  // for efficiency.
  
  const G4Axis2Placement3D& pos = *(GetPosition());
  G4Point3D p= pos.GetToPlacementCoordinates() * cp.GetPoint();

  v=forTangent*p.y()*pos.GetPX() + p.x()*pos.GetPY();
  if(GetSameSense())
    v = -v;
  
  return true;
}

