// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidCone.cc,v 1.4 2000-11-08 14:22:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidCone.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolidCone.hh"
#include "G4FPlane.hh"
#include "G4FConicalSurface.hh"
#include "G4FCylindricalSurface.hh"
#include "G4CircularCurve.hh"

G4BREPSolidCone::G4BREPSolidCone(const G4String& name,
				 const G4ThreeVector& origin,
				 const G4ThreeVector& axis,
				 const G4ThreeVector& direction,
				 G4double length,
				 G4double radius,
				 G4double large_radius)
  : G4BREPSolid(name)
{
  SurfaceVec              = new G4Surface*[3];
  G4Point3D    ArcStart1  = G4Point3D(origin + (radius*direction));
  G4Vector3D   tmpaxis(axis);
  G4Vector3D   tmporigin(origin);  
  G4Point3D    tmppoint;

  tmppoint= origin + (length*tmpaxis);
  G4Point3D origin2(tmppoint.x(), tmppoint.y(), tmppoint.z());

  tmppoint=  origin2 + (large_radius*tmpaxis);
  G4Point3D ArcStart2(tmppoint.x(), tmppoint.y(), tmppoint.z());

  G4Ray::Vcross(tmpaxis, axis, direction);
  G4ThreeVector axis2(tmpaxis.x(),tmpaxis.y(), tmpaxis.z());

  G4CurveVector CVec;
  G4CircularCurve* tmp;

  tmp = new G4CircularCurve();
  tmp->Init(G4Axis2Placement3D(direction, axis2, origin) , large_radius);
  tmp->SetBounds(ArcStart1, ArcStart1);
  CVec.insert(tmp);

  tmp = new G4CircularCurve();
  tmp->Init(G4Axis2Placement3D(direction, axis2, origin2), large_radius);
  tmp->SetBounds(ArcStart2, ArcStart2);
  CVec.insert(tmp);

  SurfaceVec[0] = new G4FConicalSurface(tmporigin, axis, 
					length, radius, large_radius);
  SurfaceVec[0]->SetBoundaries(&CVec);

  // new G4AdvancedFace("G4FConicalSurface", tmporigin, direction, 
  // axis, CVec, 1, 0,0,length, radius, large_radius);

  // Create end planes & boundaries for cone solid
  G4CurveVector CVec2;
  tmp = new G4CircularCurve();
  tmp->Init(G4Axis2Placement3D(direction, axis2, origin), radius);
  tmp->SetBounds(ArcStart1, ArcStart1);
  CVec2.insert(tmp);

  SurfaceVec[1] = new G4FPlane(tmpaxis, direction, origin2);
  //new G4AdvancedFace("G4FPlane" , origin2, direction, tmpaxis, CVec2, 1);
  SurfaceVec[1]->SetBoundaries(&CVec2);

  CVec2[0] = tmp = new G4CircularCurve();
  tmp->Init(G4Axis2Placement3D(direction, axis2, origin2), large_radius);
  tmp->SetBounds(ArcStart2, ArcStart2);

  SurfaceVec[2] = new G4FPlane(tmpaxis, direction, origin);
  //new G4AdvancedFace("G4FPlane", origin, direction, tmpaxis, CVec2, 1);  
  SurfaceVec[2]->SetBoundaries(&CVec2);

  nb_of_surfaces = 3;
  active=1;
  Initialize();
}

G4BREPSolidCone::~G4BREPSolidCone()
{
}

void G4BREPSolidCone::Initialize()
{
  // Calc bounding box for solids and surfaces
  // Convert concave planes to convex     
  ShortestDistance=1000000;
  CheckSurfaceNormals();
  if(!Box || !AxisBox)
    IsConvex();
  CalcBBoxes();
}

EInside G4BREPSolidCone::Inside(register const G4ThreeVector& Pt) const
{
  G4double dist1 = SurfaceVec[0]->HowNear(Pt);
  G4double dist2 = SurfaceVec[1]->ClosestDistanceToPoint(Pt);
  G4double dist3 = SurfaceVec[2]->ClosestDistanceToPoint(Pt);  
  if(dist1 > dist2) dist1 = dist2;
  if(dist1 > dist3) dist1 = dist3;  
  if(dist1 > 0) return kInside;
  if(dist1 < 0) return kOutside;
  return kSurface;
}

G4ThreeVector G4BREPSolidCone::SurfaceNormal(const G4ThreeVector& Pt) const
{
  G4Vector3D n =  SurfaceVec[0]->Normal(Pt);
  G4ThreeVector norm(n.x(), n.y(), n.z());
  return norm;
}

G4double G4BREPSolidCone::DistanceToIn(const G4ThreeVector& Pt) const
{
  G4double dist1 = fabs(SurfaceVec[0]->HowNear(Pt));
  G4double dist2 = fabs(SurfaceVec[1]->ClosestDistanceToPoint(Pt));
  G4double dist3 = fabs(SurfaceVec[2]->ClosestDistanceToPoint(Pt));  
  if(dist1 > dist2) dist1 = dist2;
  if(dist1 > dist3) dist1 = dist3;  
  return dist1;
 
}

G4double G4BREPSolidCone::DistanceToIn(register const G4ThreeVector& Pt, 
				       register const G4ThreeVector& V) const
{
  Reset();  
  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  //  G4double kInfinity = 10e20;
  G4Ray r(Pttmp, Vtmp);

  if(SurfaceVec[0]->Intersect( r ))
  {
    ShortestDistance = SurfaceVec[0]->GetDistance();
    return ShortestDistance;
  }
  return kInfinity; 
}

G4double G4BREPSolidCone::DistanceToOut(register const G4ThreeVector& Pt, 
					register const G4ThreeVector& V, 
					const G4bool calcNorm, 
					G4bool *validNorm, 
					G4ThreeVector *n) const
{
  if(validNorm)
    *validNorm = false;
  Reset();  

  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  //  G4double kInfinity = 10e20;

  G4Ray r(Pttmp, Vtmp);
  if(SurfaceVec[0]->Intersect( r ))
  {
    ShortestDistance = SurfaceVec[0]->GetDistance();
    return ShortestDistance;
  }
  return kInfinity; 
}

G4double G4BREPSolidCone::DistanceToOut(const G4ThreeVector& Pt) const
{
  G4double dist1 = fabs(SurfaceVec[0]->HowNear(Pt));
  G4double dist2 = fabs(SurfaceVec[1]->ClosestDistanceToPoint(Pt));
  G4double dist3 = fabs(SurfaceVec[2]->ClosestDistanceToPoint(Pt));  
  if(dist1 > dist2) dist1 = dist2;
  if(dist1 > dist3) dist1 = dist3;  
  return dist1;
}
