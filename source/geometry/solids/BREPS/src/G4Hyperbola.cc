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
// $Id: G4Hyperbola.cc,v 1.13 2010-07-07 14:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Hyperbola.cc
//
// ----------------------------------------------------------------------

#include "G4Hyperbola.hh"
#include "G4CurvePoint.hh"
#include "G4GeometryTolerance.hh"

G4Hyperbola::G4Hyperbola()
  : semiAxis(0.), semiImagAxis(0.), ratioAxisImagAxis(0.), forTangent(0.)
{
}

G4Hyperbola::~G4Hyperbola()
{
}

G4Hyperbola::G4Hyperbola(const G4Hyperbola& right)
  : G4Conic(), Focus1(right.Focus1), Focus2(right.Focus2),
    ProjFocus1(right.ProjFocus1), ProjFocus2(right.ProjFocus2),
    semiAxis(right.semiAxis), semiImagAxis(right.semiImagAxis),
    ratioAxisImagAxis(right.ratioAxisImagAxis),
    toUnitHyperbola(right.toUnitHyperbola), forTangent(right.forTangent)
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

G4Hyperbola& G4Hyperbola::operator=(const G4Hyperbola& right)
{
  if (&right == this) return *this;

  Focus1 = right.Focus1;
  Focus2 = right.Focus2;
  ProjFocus1   = right.ProjFocus1;
  ProjFocus2   = right.ProjFocus2;
  semiAxis     = right.semiAxis;
  semiImagAxis = right.semiImagAxis;
  ratioAxisImagAxis = right.ratioAxisImagAxis;
  toUnitHyperbola   = right.toUnitHyperbola;
  forTangent        = right.forTangent;
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

G4Curve* G4Hyperbola::Project(const G4Transform3D& tr)
{
  G4Exception("G4Hyperbola::Project()", "NotImplemented",
              FatalException, "Sorry, not yet implemented.");

  G4Point3D newLocation= tr*position.GetLocation();
  newLocation.setZ(0);
  G4double axisZ= (tr*position.GetPZ()).unit().z();

  if (std::abs(axisZ)<G4GeometryTolerance::GetInstance()->GetAngularTolerance()) 
  {
    return 0;
  }

  G4Vector3D newAxis(0, 0, axisZ>0? +1: -1);

  // get the parameter of an endpoint of an axis
  // (this is a point the distance of which from the center is extreme)
  G4Vector3D xPrime = tr*position.GetPX();
  xPrime.setZ(0);
  
  G4Vector3D yPrime = tr*position.GetPY();
  yPrime.setZ(0);
  
  G4Vector3D a = G4Vector3D( semiAxis*xPrime );
  G4Vector3D b = G4Vector3D( semiImagAxis*yPrime );

  G4double xval = -2*a*b/(a.mag2()+b.mag2());
  
  G4double u= (0.5*std::log((1+xval)/(1-xval)))/2;  // atanh(xval)/2

  // get the coordinate axis directions and the semiaxis lengths
  G4Vector3D sAxis= G4Vector3D( a*std::cosh(u)+b*std::sinh(u) );
  
  //!!!!!!!!!!!!
  G4Vector3D sImagAxis= G4Vector3D( a*std::cosh(u+pi/2)+b*std::sinh(u+pi/2) );
  
  //!!!!!!!!!!!!
  G4double   newSemiAxis     = sAxis.mag();
  G4double   newSemiImagAxis = sImagAxis.mag();
  G4Vector3D newRefDirection = sAxis;

  // create the new hyperbola
  G4Axis2Placement3D newPosition;
  newPosition.Init(newRefDirection, newAxis, newLocation);

  G4Hyperbola* r= new G4Hyperbola;
  r->Init(newPosition, newSemiAxis, newSemiImagAxis);

  // introduce the shift in the parametrization
  // maybe the Sign must be changed?
  r->SetPShift(u);
  
  // set the bounds when necessary
  if (IsBounded()) 
    r->SetBounds(GetPStart(), GetPEnd());
  
  return r;
}


void G4Hyperbola::InitBounded()
{
  // the bbox must include the start and endpoints as well as the
  // extreme points if they lie on the curve
  bBox.Init(GetStart(), GetEnd());

  // the parameter values
  // belonging to the points with an extreme x, y and z coordinate
  for (G4int i=0; i<3; i++) 
  {
    G4double x_i= position.GetPX()(i);
    
    if (std::abs(x_i) <= 
        G4GeometryTolerance::GetInstance()->GetAngularTolerance()) 
    {
      G4double tanhu= - (semiImagAxis*position.GetPY()(i)) / (semiAxis*x_i);
      
      if (std::abs(tanhu)<=1) 
      {
        G4double u= 0.5*std::log((1+tanhu)/(1-tanhu));  // atanh(tanhu)
	if (IsPOn(u))
	  bBox.Extend(GetPoint(u));
      }
    }
  }
}

G4bool G4Hyperbola::Tangent(G4CurvePoint& cp, G4Vector3D& v)
{
  // The tangent is computed from the 3D point representation
  // for all conics. An alternaive implementation (based on
  // the parametric point) might be worthwhile adding
  // for efficiency.
  
  const G4Axis2Placement3D& pos= *(GetPosition());
  G4Point3D p= pos.GetToPlacementCoordinates() * cp.GetPoint();
  
  v= forTangent*p.y()*pos.GetPX() + p.x()*pos.GetPY();
  
  return true;
}
