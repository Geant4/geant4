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
//
// $Id: G4VisExtent.cc 102801 2017-02-22 15:17:53Z gcosmo $
//
// 
// A.Walkden 28/11/95
// G4VisExtent.cc - to return parameters useful to the drawing window
// employed by Visualization code. 

#include "G4VisExtent.hh"

#include "G4ios.hh"

G4VisExtent::G4VisExtent (G4double xmin, G4double xmax, 
                          G4double ymin, G4double ymax, 
                          G4double zmin, G4double zmax):
  fXmin(xmin), fXmax(xmax), fYmin(ymin), fYmax(ymax), fZmin(zmin), fZmax(zmax),
  fRadiusCached(false), fCentreCached(false), fRadius(0.)
{}

G4VisExtent::G4VisExtent (const G4Point3D& centre, G4double radius):
  fRadiusCached(true), fCentreCached(true),
  fRadius(radius), fCentre(centre)
{
  // Use exscribed radius ... see comments in header file.
  G4double halfSide (radius / std::sqrt (3.));
  fXmin = centre.x () - halfSide;
  fXmax = centre.x () + halfSide;
  fYmin = centre.y () - halfSide;
  fYmax = centre.y () + halfSide;
  fZmin = centre.z () - halfSide;
  fZmax = centre.z () + halfSide;
}

G4VisExtent::~G4VisExtent () {}

const G4Point3D& G4VisExtent::GetExtentCentre () const {
  if (!fCentreCached) {
    fCentre = G4Point3D (((fXmin + fXmax) / 2.),
                         ((fYmin + fYmax) / 2.),
                         ((fZmin + fZmax) / 2.));
    fCentreCached = true;
  }
  return fCentre;
}

G4double G4VisExtent::GetExtentRadius () const {
  if (!fRadiusCached) {
    fRadius = std::sqrt (((fXmax - fXmin) * (fXmax - fXmin)) +
                         ((fYmax - fYmin) * (fYmax - fYmin)) +
                         ((fZmax - fZmin) * (fZmax - fZmin))) / 2.;
    fRadiusCached = true;
  }
  return fRadius;
}

std::ostream& operator << (std::ostream& os, const G4VisExtent& e) {
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
