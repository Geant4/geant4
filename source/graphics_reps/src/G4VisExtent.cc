// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisExtent.cc,v 1.2 1999-01-08 16:32:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// A.Walkden 28/11/95
// G4VisExtent.cc - to return parameters useful to the drawing window
// employed by Visualization code. 

#include "G4VisExtent.hh"

#include "G4ios.hh"

G4VisExtent::G4VisExtent (G4double xmin, G4double xmax, 
			  G4double ymin, G4double ymax, 
			  G4double zmin, G4double zmax) 
:fXmin(xmin), fXmax(xmax), fYmin(ymin), fYmax(ymax), fZmin(zmin), fZmax(zmax)
{}

G4VisExtent::G4VisExtent (const G4Point3D& centre, G4double radius) {
  // Use inscribed radius to define G3VisExtent so that
  // GetExtentRadius gets radius back again.  The one is the "inverse"
  // of the other, so to speak.
  G4double inscribedRadius = radius / sqrt (3.);
  fXmin = centre.x () - inscribedRadius;
  fXmax = centre.x () + inscribedRadius;
  fYmin = centre.y () - inscribedRadius;
  fYmax = centre.y () + inscribedRadius;
  fZmin = centre.z () - inscribedRadius;
  fZmax = centre.z () + inscribedRadius;
}

G4VisExtent::~G4VisExtent () {}

G4Point3D G4VisExtent::GetExtentCentre () const {
  return G4Point3D (((fXmin + fXmax) / 2.),
		       ((fYmin + fYmax) / 2.),
		       ((fZmin + fZmax) / 2.));
}

G4double G4VisExtent::GetExtentRadius () const {
  return sqrt (((fXmax - fXmin) * (fXmax - fXmin)) +
	       ((fYmax - fYmin) * (fYmax - fYmin)) +
	       ((fZmax - fZmin) * (fZmax - fZmin))) / 2.;
}
 
ostream& operator << (ostream& os, const G4VisExtent& e) {
  os << "G4VisExtent (bounding box):";
  os << "\n  X limits: " << e.fXmin << ' ' << e.fXmax;
  os << "\n  Y limits: " << e.fYmin << ' ' << e.fYmax;
  os << "\n  Z limits: " << e.fZmin << ' ' << e.fZmax;
  return os;
}

G4bool operator != (const G4VisExtent& e1, const G4VisExtent& e2) {
  return ((e1.fXmin != e2.fXmin) ||
	  (e1.fXmax != e2.fXmax) ||
	  (e1.fYmin != e2.fYmin) ||
	  (e1.fYmax != e2.fYmax) ||
	  (e1.fZmin != e2.fZmin) ||
	  (e1.fZmax != e2.fZmax));
}
