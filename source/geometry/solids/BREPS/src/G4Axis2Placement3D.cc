// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Axis2Placement3D.cc,v 1.4 2000-11-10 17:41:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Axis2Placement3D.cc
//
// ----------------------------------------------------------------------

#include "G4Axis2Placement3D.hh"

//G4Axis2Placement3D
G4Axis2Placement3D::G4Axis2Placement3D(){}
G4Axis2Placement3D::~G4Axis2Placement3D(){}

// copy constructor (used in STEPinterface module)
//
G4Axis2Placement3D::G4Axis2Placement3D(const G4Axis2Placement3D& place)
  : refDirection(place.refDirection), 
    axis(place.axis), location(place.location),
    pX(place.pX), pY(place.pY), pZ(place.pZ),
    toPlacementCoordinates(place.toPlacementCoordinates),
    fromPlacementCoordinates(place.fromPlacementCoordinates)
{
}

// assignment operator
//
G4Axis2Placement3D&
G4Axis2Placement3D::operator=(const G4Axis2Placement3D& place)
{
  if (&place == this) return *this;
  
  refDirection             = place.refDirection; 
  axis                     = place.axis;
  location                 = place.location;
  pX                       = place.pX;
  pY                       = place.pY;
  pZ                       = place.pZ;
  toPlacementCoordinates   = place.toPlacementCoordinates;
  fromPlacementCoordinates = place.fromPlacementCoordinates; 

  return *this;
}

/* everything below here is commented-out ...

G4Axis2Placement3D::G4Axis2Placement3D(const G4ThreeVec Dir, 
                                       const G4ThreeVec Axis, 
				       const G4Point3d Pt    )
{
  dir=Dir;
  axis=Axis;
  srf_point=Pt;
  ComputeNormal();
  G4Point3d Pt2 = Pt+Dir;
  G4Point3d Pt3 = Pt+Axis;  
  G4Ray::CalcPlane3Pts(Pl, Pt, Pt2, Pt3);  
}

G4Axis2Placement3D::G4Axis2Placement3D(const G4ThreeVec Dir,
                                       const G4ThreeVec Axis,
				       const G4Point3d Pt1,
				       const G4Point3d Pt2,
				       const G4Point3d Pt3)
{
  dir=Dir;
  axis=Axis;
  srf_point=Pt1;
  ComputeNormal();
  G4Ray::CalcPlane3Pts(Pl, Pt1, Pt2, Pt3);  
}

void
G4Axis2Placement3D::ProjectPlacement(const G4Plane& Pl1,
                                     const G4Plane& Pl2)
{
  Project(ProjectedDir, dir, Pl1, Pl2);
  Project(ProjectedAxis, axis, Pl1, Pl2);
  Project(ProjectedSrfPoint, srf_point, Pl1, Pl2);
  Project(ProjectedNormal, Normal, Pl1, Pl2);
}

void
G4Axis2Placement3D::ComputeNormal()
{

  if(dir == axis)
    Normal = dir;
  else
    {
      Normal.X(dir.Y()*axis.Z() - dir.Z()*axis.Y());
      Normal.Y(dir.X()*axis.Z()- dir.Z()*axis.X());
      Normal.Z(dir.X()*axis.Y() - dir.Y()*axis.X());
    }
}


G4Point3d
G4Axis2Placement3D::EvaluateIntersection(register const G4Ray& rray)
{

// s is solution, line is p + tq, n is G4Plane Normal, r is point on G4Plane 
// all parameters are pointers to arrays of three elements 

    register G4double a, b, t;
    register const G4ThreeVec& RayDir   = rray.GetDir();
    register const G4Point3d& RayStart = rray.GetStart();
    G4double dirx =  RayDir.X();
    G4double diry =  RayDir.Y();
    G4double dirz =  RayDir.Z();
    b = Normal.X() * dirx + Normal.Y() * diry + Normal.Z() * dirz;

    if (fabs(b) < 0.001)//== 0.0)
       // or some better test involving a small positive e     
    {
//    G4cout << "\nLine is parallel to G4Plane.No Hit.";
      G4Point3d hit_point( kInfinity, kInfinity, kInfinity);
      return hit_point;
    }
    G4double startx =  RayStart.X();
    G4double starty =  RayStart.Y();
    G4double startz =  RayStart.Z();    
    
    a = Normal.X() * (srf_point.X() - startx)
      + Normal.Y() * (srf_point.Y() - starty)
      + Normal.Z() * (srf_point.Z() - startz);

    t = a/b;

    // substitute t into line equation
    // to calculate final solution     
    G4Point3d hit_point(startx + t * dirx,starty
                               + t * diry,startz
			       + t * dirz);    
    
//  G4cout << "\nPLANE HIT POINT :" << hit_point.X()
//         << "  " << hit_point.Y() << "  " << hit_point.Z();
    return hit_point;
}
*/
