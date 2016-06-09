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
// G4CircularCurve.cc
//
// ----------------------------------------------------------------------

#include "G4CircularCurve.hh"
#include "G4PhysicalConstants.hh"
#include "G4Ellipse.hh"

// G4CircularCurve
G4CircularCurve::G4CircularCurve() : radius(0.) {}
G4CircularCurve::~G4CircularCurve() {}

G4CircularCurve::G4CircularCurve(const G4CircularCurve& right)
  : G4Conic(), radius(right.radius)
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

G4CircularCurve&
G4CircularCurve::operator=(const G4CircularCurve& right)
{
  if (&right == this) return *this;

  radius    = right.radius;
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

  return *this;
}

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
    G4double u = std::atan2(position.GetPY()(i), position.GetPX()(i));
    
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

G4double G4CircularCurve::GetPMax() const
{
  return twopi;
}

G4Point3D G4CircularCurve::GetPoint(G4double param) const
{
  return G4Point3D( position.GetLocation()+radius*
    ( std::cos(param)*position.GetPX() + std::sin(param)*position.GetPY() ) );
}

G4double G4CircularCurve::GetPPoint(const G4Point3D& pt) const
{
  G4Point3D ptLocal= position.GetToPlacementCoordinates()*pt;
  G4double angle= std::atan2(ptLocal.y(), ptLocal.x());
  return (angle<0)? angle+twopi: angle;
}

////////////////////////////////////////////////////////////////////////////

/*
#include "G4CurveRayIntersection.hh"

void G4CircularCurve::IntersectRay2D(const G4Ray&,
				     G4CurveRayIntersection&)
{
  G4Exception("G4CircularCurve::IntersectRay2D()",
              "NotApplicable", FatalException,
              "G4CircularCurve is always 3D!");
}
*/

G4int G4CircularCurve::IntersectRay2D(const G4Ray&)
{
  G4Exception("G4CircularCurve::IntersectRay2D()", "GeomSolids0002",
              FatalException,"G4CircularCurve is always 3D!");
  return 0;
}
