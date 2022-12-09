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
//
// 
// A.Walkden 28/11/95
// G4VisExtent.cc - to return parameters useful to the drawing window
// employed by Visualization code. 

#include "G4VisExtent.hh"

#include "G4ThreeVector.hh"

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

G4VisExtent::~G4VisExtent () = default;

const G4VisExtent& G4VisExtent::GetNullExtent () {
  static const G4VisExtent nullExtent = G4VisExtent();
  return nullExtent;
}

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

G4VisExtent& G4VisExtent::Transform (const G4Transform3D& transform)
{
  const auto& rotation = transform.getRotation();
  const auto& translation = transform.getTranslation();

  G4ThreeVector nnn(fXmin,fYmin,fZmin);
  G4ThreeVector nnx(fXmin,fYmin,fZmax);
  G4ThreeVector nxn(fXmin,fYmax,fZmin);
  G4ThreeVector nxx(fXmin,fYmax,fZmax);
  G4ThreeVector xnn(fXmax,fYmin,fZmin);
  G4ThreeVector xnx(fXmax,fYmin,fZmax);
  G4ThreeVector xxn(fXmax,fYmax,fZmin);
  G4ThreeVector xxx(fXmax,fYmax,fZmax);

  nnn.transform(rotation); nnn += translation;
  nnx.transform(rotation); nnx += translation;
  nxn.transform(rotation); nxn += translation;
  nxx.transform(rotation); nxx += translation;
  xnn.transform(rotation); xnn += translation;
  xnx.transform(rotation); xnx += translation;
  xxn.transform(rotation); xxn += translation;
  xxx.transform(rotation); xxx += translation;

  fXmin = DBL_MAX;
  fYmin = DBL_MAX;
  fZmin = DBL_MAX;
  fXmax = -DBL_MAX;
  fYmax = -DBL_MAX;
  fZmax = -DBL_MAX;
  for (const auto& corner: {nnn,nnx,nxn,nxx,xnn,xnx,xxn,xxx}) {
    if (fXmin > corner.getX()) fXmin = corner.getX();
    if (fYmin > corner.getY()) fYmin = corner.getY();
    if (fZmin > corner.getZ()) fZmin = corner.getZ();
    if (fXmax < corner.getX()) fXmax = corner.getX();
    if (fYmax < corner.getY()) fYmax = corner.getY();
    if (fZmax < corner.getZ()) fZmax = corner.getZ();
  }

  return *this;
}
