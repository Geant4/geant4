// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidCylinder.cc,v 1.4 2000-11-08 14:22:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidCylinder.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolidCylinder.hh"
#include "G4CircularCurve.hh"
#include "G4FPlane.hh"
#include "G4FCylindricalSurface.hh"

G4BREPSolidCylinder::G4BREPSolidCylinder(const G4String& name,
					 const G4ThreeVector& origin,
					 const G4ThreeVector& axis,
					 const G4ThreeVector& direction,
					 G4double radius,
					 G4double length)
  : G4BREPSolid(name)
{
  SurfaceVec = new G4Surface*[3];
  G4CurveVector cv;
  G4CircularCurve* tmp;



  // Creation of the cylindrical surface
  SurfaceVec[0] = new G4FCylindricalSurface(origin, axis, radius , length);
  //SurfaceVec[0]->SetBoundaries(&cv);
  //cv.clear();


  // Creation of the first circular surface, which origin is origin
  G4Point3D  ArcStart1 = G4Point3D( origin + ( radius*direction ) );
  G4Vector3D axis1     = G4Vector3D( axis.cross( direction ) );

  tmp = new G4CircularCurve;
  tmp->Init( G4Axis2Placement3D(direction, axis1, origin), radius );
  tmp->SetBounds(ArcStart1, ArcStart1);
  cv.insert(tmp);

  SurfaceVec[1] = new G4FPlane(direction, axis1, origin);
  SurfaceVec[1]->SetBoundaries(&cv);
  cv.clear();
  

  // Creation of the second circular surface
  G4Point3D  origin2   = G4Point3D( origin  + ( length*axis ) );  
  G4Point3D  ArcStart2 = G4Point3D( origin2 + ( radius*direction ) );
  G4Vector3D axis2     = axis1;

  tmp = new G4CircularCurve;
  tmp->Init( G4Axis2Placement3D(direction, axis2, origin2), radius);
  tmp->SetBounds(ArcStart2, ArcStart2);
  cv.insert(tmp);

  SurfaceVec[2] = new G4FPlane(direction, axis2, origin2);
  SurfaceVec[2]->SetBoundaries(&cv);
  cv.clear();


  nb_of_surfaces = 3;
  active=1;
  Initialize();
}

G4BREPSolidCylinder::~G4BREPSolidCylinder()
{
}
