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
//////////////////////////////////////////////////////////////////////////
// $Id: CurveTestFunction.hh,v 1.9 2007-05-18 10:31:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//////////////////////////////////////////////////////////////////////////
//
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
#include "Randomize.hh"



// Printing debug info
std::ostream& operator<<(std::ostream& os, const G4Axis2Placement3D& p)
{
  os << "(" << p.GetLocation()
     << ", " << p.GetAxis()
     << ", " << p.GetRefDirection() << ")";
  return os;
}


std::ostream& operator<<(std::ostream& os, const G4BoundingBox3D& b)
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
  
  G4double r= G4RandFlat::shoot(min, max);
  
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
      u2-= std::floor((u2-u)/pMax+0.5)*pMax;
    
    G4bool error= std::abs(u-u2) > 0.0001;
    
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

  for (i=0; i<1000; i++) 
  {
    G4double  u  = GetRandomP(c);
    G4Point3D p1 = tr*c->GetPoint(u);
    p1.setZ(0);

    G4Point3D p2 = projC->GetPoint(u);
    G4bool error = false;
 
    if (std::abs(p2.z())>0.0001) 
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
