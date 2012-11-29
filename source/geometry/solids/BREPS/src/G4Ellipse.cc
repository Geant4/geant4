//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Ellipse.cc
//
// ----------------------------------------------------------------------

#include "G4Ellipse.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"

G4Ellipse::G4Ellipse()
  : semiAxis1(0.), semiAxis2(0.), ratioAxis2Axis1(0.), forTangent(.0)
{
}

G4Ellipse::~G4Ellipse()
{
}

G4Ellipse::G4Ellipse(const G4Ellipse& right)
  : G4Conic(), semiAxis1(right.semiAxis1), semiAxis2(right.semiAxis2),
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

  if (std::abs(axisZ)<G4GeometryTolerance::GetInstance()->GetAngularTolerance()) 
    { return 0; }
  
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
  else if ((std::abs(abmag) < perMillion) && (std::abs(prod) < perMillion))
    u = 0.;
  else
    u = std::atan2(prod,abmag) / 2;

  // get the coordinate axis directions and the semiaxis lengths
  G4Vector3D sAxis1          = G4Vector3D( a*std::cos(u)+b*std::sin(u) );
  G4Vector3D sAxis2          = G4Vector3D( a*std::cos(u+pi/2)+b*std::sin(u+pi/2) );
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
  // G4double maxEntent      = std::max(semiAxis1, semiAxis2);
  // G4Vector3D halfExtent(maxEntent, maxEntent, maxEntent);
  // bBox.Init(center+halfExtent, center-halfExtent);

  // the bbox must include the start and endpoints as well as the
  // extreme points if they lie on the curve
  bBox.Init(GetStart(), GetEnd());

  // the parameter values
  // belonging to the points with an extreme x, y and z coordinate
  for (G4int i=0; i<3; i++) 
  {
    G4double u= std::atan2(position.GetPY()(i)*semiAxis2,
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

