// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Line.cc,v 1.2 1999-01-14 16:06:46 broglia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4Line.hh"

G4Line::G4Line (){}
G4Line::~G4Line (){}

G4Curve* G4Line::Project(const G4Transform3D& tr)
{
  G4Vector3D newDir= tr*dir;
  
  if (abs(newDir.x())+abs(newDir.y()) < kCarTolerance) 
    return 0;
  
  G4Point3D newPnt= tr*pnt;
  newDir.setZ(0);
  newPnt.setZ(0);
  
  G4Line* r= new G4Line();

  // L. Broglia : terrible mistake !!!!
  //r->Init(newDir, newPnt);
  r->Init(newPnt, newDir);

  r->SetBounds(GetPStart(), GetPEnd());
  
  return r;
}

////////////////////////////////////////////////////////////////////////////

G4bool G4Line::Tangent(G4CurvePoint& cp, G4Vector3D& vec)
{
  if(GetSameSense())
    vec = -dir;
  else
    vec = dir;

  return true;
}

