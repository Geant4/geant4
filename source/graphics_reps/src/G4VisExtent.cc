//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VisExtent.cc,v 1.6 2001-07-24 21:39:50 johna Exp $
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
  // Use exscribed radius ... see comments in header file.
  G4double halfSide (radius / sqrt (3.));
  fXmin = centre.x () - halfSide;
  fXmax = centre.x () + halfSide;
  fYmin = centre.y () - halfSide;
  fYmax = centre.y () + halfSide;
  fZmin = centre.z () - halfSide;
  fZmax = centre.z () + halfSide;
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
 
G4std::ostream& operator << (G4std::ostream& os, const G4VisExtent& e) {
  os << "G4VisExtent (bounding box):";
  os << "\n  X limits: " << e.fXmin << ' ' << e.fXmax;
  os << "\n  Y limits: " << e.fYmin << ' ' << e.fYmax;
  os << "\n  Z limits: " << e.fZmin << ' ' << e.fZmax;
  return os;
}

G4bool G4VisExtent::operator != (const G4VisExtent& e) const {
  return ((fXmin != e.fXmin) ||
	  (fXmax != e.fXmax) ||
	  (fYmin != e.fYmin) ||
	  (fYmax != e.fYmax) ||
	  (fZmin != e.fZmin) ||
	  (fZmax != e.fZmax));
}
