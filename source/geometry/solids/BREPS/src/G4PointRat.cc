// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PointRat.cc,v 1.3 2000-01-21 13:47:51 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Modif 8 oct 98 : A.Floquet
//      G4PointRat datas are made of
// 	     . a point 3D
//	     . a additional value : the scale factor which is set to 1 by default
//

#include "G4PointRat.hh"

G4PointRat::G4PointRat():pt3d(){s=1;}

G4PointRat::G4PointRat(const G4Point3D& tmp):pt3d(tmp){s=1;}

G4PointRat::~G4PointRat(){}

void G4PointRat::operator=(const G4PointRat& a)
{   pt3d.setX(a.x());
    pt3d.setY(a.y());
    pt3d.setZ(a.z());
    s=a.w();
}

void G4PointRat::operator=(const G4Point3D& a)
{   pt3d = a;
    s=1;
}
