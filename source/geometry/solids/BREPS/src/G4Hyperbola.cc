// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Hyperbola.cc,v 1.4 2000-08-28 15:00:38 gcosmo Exp $
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

G4Hyperbola::G4Hyperbola()
{
}

G4Hyperbola::~G4Hyperbola()
{
}

G4Curve* G4Hyperbola::Project(const G4Transform3D& tr)
{
  G4Exception("G4Hyperbola::Project");

  G4Point3D newLocation= tr*position.GetLocation();
  newLocation.setZ(0);
  G4double axisZ= (tr*position.GetPZ()).unit().z();

  if (abs(axisZ)<kAngTolerance) 
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
  
  G4Vector3D a = semiAxis*xPrime;
  G4Vector3D b = semiImagAxis*yPrime;

  G4double xval = -2*a*b/(a.mag2()+b.mag2());
  
#ifdef WIN32
  G4double u= (0.5*log((1+xval)/(1-xval)))/2;
#else
  G4double u= atanh( xval ) / 2;
#endif

  // get the coordinate axis directions and the semiaxis lengths
  G4Vector3D sAxis= a*cosh(u)+b*sinh(u);
  
  //!!!!!!!!!!!!
  G4Vector3D sImagAxis= a*cosh(u+pi/2)+b*sinh(u+pi/2);
  
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
    
    if (abs(x_i) <= kAngTolerance) 
    {
      G4double tanhu= - (semiImagAxis*position.GetPY()(i)) / (semiAxis*x_i);
      
      if (abs(tanhu)<=1) 
      {
        #ifdef WIN32
          G4double u= 0.5*log((1+tanhu)/(1-tanhu));
        #else
	  G4double u= atanh(tanhu);
        #endif
	
	  if (IsPOn(u)) bBox.Extend(GetPoint(u));
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
