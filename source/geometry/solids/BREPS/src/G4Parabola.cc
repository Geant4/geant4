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
// $Id: G4Parabola.cc,v 1.10 2010-07-07 14:45:31 gcosmo Exp $
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
#include "G4GeometryTolerance.hh"

G4Parabola::G4Parabola() : focalDist(0.) {}
G4Parabola::~G4Parabola(){}

G4Parabola::G4Parabola(const G4Parabola& right)
  : G4Conic(), focalDist(right.focalDist), F(right.F), L0(right.L0)
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

G4Parabola& G4Parabola::operator=(const G4Parabola& right)
{
  if (&right == this) return *this;

  F  = right.F;
  L0 = right.L0;
  focalDist = right.focalDist;
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

G4Curve* G4Parabola::Project(const G4Transform3D& tr)
{
  G4double axisZ= (tr*position.GetPZ()).unit().z();

  if (std::abs(axisZ)<G4GeometryTolerance::GetInstance()->GetAngularTolerance()) 
    { return 0; }
  
  
  G4Vector3D newAxis(0, 0, axisZ>0? +1: -1);

  G4Vector3D xPrime= tr*position.GetPX();
  xPrime.setZ(0);
  G4Vector3D yPrime= tr*position.GetPY();
  yPrime.setZ(0);
  G4double u= -(xPrime*yPrime)/xPrime.mag2();

  G4Point3D newLocation= G4Point3D( tr*position.GetLocation()+
                                    focalDist*(u*u*xPrime+2*u*yPrime) );
  newLocation.setZ(0);
  G4Vector3D newRefDirection= xPrime;
  G4double newFocalDist= (focalDist*((2*u+1)*xPrime+2*yPrime)).mag()/std::sqrt(5.);
			  
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
    
    if (std::abs(x_i) <= 
        G4GeometryTolerance::GetInstance()->GetAngularTolerance()) 
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
