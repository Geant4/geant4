// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidBox.cc,v 1.3 2000-08-28 08:57:55 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidBox.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolidBox.hh"
#include "G4FPlane.hh"
#include "G4Point3DVector.hh"

G4BREPSolidBox::G4BREPSolidBox(const G4String& name,
			       const G4Point3D& Pt1,
			       const G4Point3D& Pt2,
			       const G4Point3D& Pt3,
			       const G4Point3D& Pt4,
			       const G4Point3D& Pt5,
			       const G4Point3D& Pt6,
			       const G4Point3D& Pt7,
			       const G4Point3D& Pt8): G4BREPSolid(name)
{
  nb_of_surfaces=6;
  active=1;PlaneSolid=1;
  SurfaceVec = new G4Surface*[6];
  G4Point3DVector PVec(4);

  PVec[0] = Pt1;
  PVec[1] = Pt2;
  PVec[2] = Pt3;
  PVec[3] = Pt4;  
  SurfaceVec[0] = new G4FPlane(&PVec);

  PVec[2] = Pt6;
  PVec[3] = Pt5;  
  SurfaceVec[1] = new G4FPlane(&PVec);

  PVec[0] = Pt2;
  PVec[1] = Pt6;
  PVec[2] = Pt7;
  PVec[3] = Pt3;  
  SurfaceVec[2] = new G4FPlane(&PVec);

  PVec[0] = Pt3;
  PVec[1] = Pt7;
  PVec[2] = Pt8;
  PVec[3] = Pt4;  
  SurfaceVec[3] = new G4FPlane(&PVec);

  PVec[0] = Pt1;
  PVec[1] = Pt5;
  PVec[2] = Pt8;
  PVec[3] = Pt4;  
  SurfaceVec[4] = new G4FPlane(&PVec);

  PVec[0] = Pt5;
  PVec[1] = Pt6;
  PVec[2] = Pt7;
  PVec[3] = Pt8;  
  SurfaceVec[5] = new G4FPlane(&PVec);

  Initialize();
  
}

EInside G4BREPSolidBox::Inside(register const G4ThreeVector& Pt) const
{
  G4Point3D Point(Pt);

  // Get the bounding box extent
  G4Point3D min = bbox->GetBoxMin();
  min = min + -(0.5*kCarTolerance);

  G4Point3D max = bbox->GetBoxMax();
  max = max + (0.5*kCarTolerance);

  if( (Point.x() < min.x() || Point.x() > max.x()) ||
      (Point.y() < min.y() || Point.y() > max.y()) ||
      (Point.z() < min.z() || Point.z() > max.z())    )
    return kOutside;

  if( (Point.x() > min.x() && Point.x() < max.x())&&
      (Point.y() > min.y() && Point.y() < max.y())&&
      (Point.z() > min.z() && Point.z() < max.z())   )
    return kInside;
  
  return kSurface;
}
