//////////////////////////////////////////////////////////////////////////
// $Id: CurveTestFunction.hh,v 1.3 2000-08-28 08:58:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//////////////////////////////////////////////////////////////////////////
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// L. Broglia
// October 26, 1998
//
// Functions wich tested the basic functionnality of the different curves :
//   - line
//   - circular curve
//   - ellipse
//   - parabola
//   - hyperbola
//
// Test the functionality of Axis2Placement3D and Project 
//
// These functions are utilized into CurveTest.cc
//


#include "G4Axis2Placement3D.hh"
#include "G4Curve.hh"
#include "G4Line.hh"
#include "G4CircularCurve.hh"
#include "G4Ellipse.hh"
#include "G4Parabola.hh"
#include "G4Hyperbola.hh"
#include "G4Transform3D.hh"
#include "CLHEP/Random/RandFlat.h"



// Printing debug info
G4std::ostream& operator<<(G4std::ostream& os, const G4Axis2Placement3D& p)
{
  os << "(" << p.GetLocation()
     << ", " << p.GetAxis()
     << ", " << p.GetRefDirection() << ")";
  return os;
}


G4std::ostream& operator<<(G4std::ostream& os, const G4BoundingBox3D& b)
{
  os << "(" << b.GetBoxMin()
     << ", " << b.GetBoxMax() << ")";
  return os;
}


// Get a random parametric point
G4double GetRandomP(G4Curve* c)
{
  G4double pMax = c->GetPMax();
  G4double min  = c->GetPStart();
  G4double max  = c->GetPEnd();

  if ( max<=min && pMax>0) 
    max+= pMax;
  
  G4double r= RandFlat::shoot(min, max);
  
  if (pMax>0 && r>=pMax) 
    r-= pMax;
  
  return r;
}


// Test basic functionality common to curves
void TestCurve(G4Curve* c)
{
  // bounds
  G4cout << "  GetStartPoint  : " << c->GetStart()  << G4endl;
  G4cout << "  GetEndPoint    : " << c->GetEnd()    << G4endl;
  G4cout << "  GetPStartPoint : " << c->GetPStart() << G4endl;
  G4cout << "  GetPEndPoint   : " << c->GetPEnd()   << G4endl;
  G4cout << "  GetPMax        : " << c->GetPMax()   << G4endl;
  G4cout << "  u -> GetPoint() -> p -> GetPPoint() -> u2" << G4endl;

  G4int i;
  for (i=0; i<1000; i++) 
  {
    G4double  u    = GetRandomP(c);
    G4Point3D p    = c->GetPoint(u);
    G4double  u2   = c->GetPPoint(p);
    G4double  pMax = c->GetPMax();

    if (pMax>0) 
      u2-= floor((u2-u)/pMax+0.5)*pMax;
    
    G4bool error= abs(u-u2) > 0.0001;
    
    if (i<5 || error) 
    {
      G4cout << " u : " << u << " p: " << p
	     << " u2: " << u2 << G4endl;
      
      if (error) 
      {
	G4cout << "    Error: u != u2 !" << G4endl;
	exit(EXIT_FAILURE);
      }
    }

    if (!c->IsPOn(u)) 
    {
      G4cout << "    Error: u is not on the curve!" << G4endl;
      exit(EXIT_FAILURE);
    }
  }

  const G4BoundingBox3D& bBox= *(c->BBox());
  G4cout << "  BBox: " << bBox << G4endl;
  G4BoundingBox3D bBox2(c->GetStart(), c->GetEnd());
  
  for (i=0; i<1000; i++) 
    bBox2.Extend(c->GetPoint(GetRandomP(c)));
 
  G4cout << "  BBox from random points: " << bBox2 << G4endl;
  G4cout << "  empty space around it when put in BBox:" << G4endl;
  
  G4Vector3D max= bBox.GetBoxMax()-bBox2.GetBoxMax();
  G4Vector3D min= bBox2.GetBoxMin()-bBox.GetBoxMin();
  
  G4cout << " min x: " << min.x() << " max x: " << max.x()
	 << "\n min y: " << min.y() << " max y: " << max.y()
	 << "\n min z: " << min.z() << " max z: " << max.z() << G4endl;
}


// Test G4Axis2Placement3D
void TestPlacement(G4Axis2Placement3D* p)
{
  G4cout << "G4Axis2Placement3D " << p->GetLocation() << " " << p->GetAxis()
	 << " " << p->GetRefDirection() << G4endl;

  G4cout << "  coordinate vectors:" << G4endl;

  G4cout << " p[1]: " << p->GetPX() 
	 << "\n p[2]: " << p->GetPY()
	 << "\n p[3]: " << p->GetPZ() << G4endl;

  G4cout << " p[1].p[2]: " << p->GetPX()*p->GetPY()
	 << "\n p[3]-(p[1]xp[2]): "<< p->GetPZ()-(p->GetPX().cross(p->GetPY()))
	 << G4endl;

  G4cout << "  coordinate vectors computed from the transformation:" << G4endl;

  G4cout << " p[1]: " << p->GetFromPlacementCoordinates()*G4Point3D(1, 0, 0)
	 << "\n p[2]: " << p->GetFromPlacementCoordinates()*G4Point3D(0, 1, 0)
	 << "\n p[3]: " << p->GetFromPlacementCoordinates()*G4Point3D(0, 0, 1)
	 << G4endl;
  
  G4cout << "  coordinate vectors in placement coordinate system:";

  G4cout << "\n   " << p->GetToPlacementCoordinates()*p->GetPX()
	 << "\n   " << p->GetToPlacementCoordinates()*p->GetPY()
	 << "\n   " << p->GetToPlacementCoordinates()*p->GetPZ() << G4endl;
}


// Test the project function
void TestProject(G4Curve* c, const G4Transform3D& tr)
{
  G4int i;
  G4Curve*   projC = c->Project(tr);
  G4Ellipse* lof   = (G4Ellipse*)projC;

  for (i=0; i<1000; i++) 
  {
    G4double  u  = GetRandomP(c);
    G4Point3D p1 = tr*c->GetPoint(u);
    p1.setZ(0);

    G4Point3D p2 = projC->GetPoint(u);
    G4bool error = false;
 
    if (abs(p2.z())>0.0001) 
    {
      error= true;
      G4cout << "Error: the projection is not in the xy plane!" << G4endl;
    }

    if ((p1-p2).mag()>0.001) 
    {
      error= true;
      G4cout << "Error: the two projected points should be the same!" << G4endl;
    }

    if (error) 
    {
      G4cout << "u: " << u << G4endl;
     
      G4cout << "point on original curve:\n" << c->GetPoint(u)
	     << "\n projected: " << tr*c->GetPoint(u) << G4endl;
      
      G4cout << "point on projected curve: " << projC->GetPoint(u) << G4endl;
      
      u= -u;
      G4cout << "u: " << u << G4endl;
      
      G4cout << "point on original curve:\n" << c->GetPoint(u)
	     << "\n projected: " << tr*c->GetPoint(u) << G4endl;
     
      G4cout << "point on projected curve:\n" << projC->GetPoint(u) << G4endl;
      
      u= -u; u+= pi;
      G4cout << "u: " << u << G4endl;
      
      G4cout << "point on original curve:\n" << c->GetPoint(u)
	     << "\n projected: " << tr*c->GetPoint(u) << G4endl;
      
      G4cout << "\npoint on projected curve: " << projC->GetPoint(u) << G4endl;
      
      u= -u;
      G4cout << "u: " << u << G4endl;
      
      G4cout << "point on original curve:\n" << c->GetPoint(u)
	     << "\n projected: " << tr*c->GetPoint(u) << G4endl;
      
      G4cout << "point on projected curve: " << projC->GetPoint(u) << G4endl;
      
      exit(EXIT_FAILURE);
    }
  }
}
