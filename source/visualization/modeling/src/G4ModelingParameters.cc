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
// $Id: G4ModelingParameters.cc,v 1.10 2006/06/29 21:32:48 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// John Allison  31st December 1997.
// Parameters associated with the modeling of GEANT4 objects.

#include "G4ModelingParameters.hh"

#include "G4ios.hh"
#include "G4VisAttributes.hh"
#include "G4ExceptionSeverity.hh"

G4ModelingParameters::G4ModelingParameters ():
  fpDefaultVisAttributes (0),
  fDrawingStyle          (wf),
  fRepStyle              (polyhedron),
  fCulling               (false),
  fCullInvisible         (false),
  fDensityCulling        (false),
  fVisibleDensity        (0.01 * g / cm3),
  fCullCovered           (false),
  fNoOfSides             (24),
  fViewGeom              (true),
  fViewHits              (true),
  fViewDigis             (true)
{}

G4ModelingParameters::G4ModelingParameters
(const G4VisAttributes* pDefaultVisAttributes,
 G4ModelingParameters::RepStyle repStyle,
 G4bool isCulling,
 G4bool isCullingInvisible,
 G4bool isDensityCulling,
 G4double visibleDensity,
 G4bool isCullingCovered,
 G4int noOfSides
 ):
  fpDefaultVisAttributes (pDefaultVisAttributes),
  fDrawingStyle   (wf),
  fRepStyle       (repStyle),
  fCulling        (isCulling),
  fCullInvisible  (isCullingInvisible),
  fDensityCulling (isDensityCulling),
  fVisibleDensity (visibleDensity),
  fCullCovered    (isCullingCovered),
  fNoOfSides      (noOfSides),
  fViewGeom       (true),
  fViewHits       (true),
  fViewDigis      (true)
{}

G4ModelingParameters::G4ModelingParameters
(const G4VisAttributes* pDefaultVisAttributes,
 G4ModelingParameters::RepStyle repStyle,
 G4bool isCulling,
 G4bool isCullingInvisible,
 G4bool isDensityCulling,
 G4double visibleDensity,
 G4bool isCullingCovered,
 G4int noOfSides,
 G4bool isViewGeom,
 G4bool isViewHits,
 G4bool isViewDigis
 ):
  fpDefaultVisAttributes (pDefaultVisAttributes),
  fDrawingStyle   (wf),
  fRepStyle       (repStyle),
  fCulling        (isCulling),
  fCullInvisible  (isCullingInvisible),
  fDensityCulling (isDensityCulling),
  fVisibleDensity (visibleDensity),
  fCullCovered    (isCullingCovered),
  fNoOfSides      (noOfSides),
  fViewGeom       (isViewGeom),
  fViewHits       (isViewHits),
  fViewDigis      (isViewDigis)
{}

G4ModelingParameters::G4ModelingParameters
(const G4VisAttributes* pDefaultVisAttributes,
 G4ModelingParameters::DrawingStyle drawingStyle,
 G4ModelingParameters::RepStyle repStyle,
 G4bool isCulling,
 G4bool isCullingInvisible,
 G4bool isDensityCulling,
 G4double visibleDensity,
 G4bool isCullingCovered,
 G4int noOfSides
 ):
  fpDefaultVisAttributes (pDefaultVisAttributes),
  fDrawingStyle   (drawingStyle),
  fRepStyle       (repStyle),
  fCulling        (isCulling),
  fCullInvisible  (isCullingInvisible),
  fDensityCulling (isDensityCulling),
  fVisibleDensity (visibleDensity),
  fCullCovered    (isCullingCovered),
  fNoOfSides      (noOfSides),
  fViewGeom       (true),
  fViewHits       (true),
  fViewDigis      (true)
{}

G4ModelingParameters::G4ModelingParameters
(const G4VisAttributes* pDefaultVisAttributes,
 G4ModelingParameters::DrawingStyle drawingStyle,
 G4ModelingParameters::RepStyle repStyle,
 G4bool isCulling,
 G4bool isCullingInvisible,
 G4bool isDensityCulling,
 G4double visibleDensity,
 G4bool isCullingCovered,
 G4int noOfSides,
 G4bool isViewGeom,
 G4bool isViewHits,
 G4bool isViewDigis
 ):
  fpDefaultVisAttributes (pDefaultVisAttributes),
  fDrawingStyle   (drawingStyle),
  fRepStyle       (repStyle),
  fCulling        (isCulling),
  fCullInvisible  (isCullingInvisible),
  fDensityCulling (isDensityCulling),
  fVisibleDensity (visibleDensity),
  fCullCovered    (isCullingCovered),
  fNoOfSides      (noOfSides),
  fViewGeom       (isViewGeom),
  fViewHits       (isViewHits),
  fViewDigis      (isViewDigis)
{}

G4ModelingParameters::~G4ModelingParameters () {}

void G4ModelingParameters::SetVisibleDensity (G4double visibleDensity) {
  const G4double reasonableMaximum = 10.0 * g / cm3;
  if (visibleDensity < 0) {
    G4cout << "G4ModelingParameters::SetVisibleDensity: attempt to set negative "
      "density - ignored." << G4endl;
  }
  else {
    if (fVisibleDensity > reasonableMaximum) {
      G4cout << "G4ModelingParameters::SetVisibleDensity: density > "
	   << reasonableMaximum
	   << " g / cm3 - did you mean this?"
	   << G4endl;
    }
    fVisibleDensity = visibleDensity;
  }
}

G4int G4ModelingParameters::SetNoOfSides (G4int nSides) {
  const G4int  nSidesMin = 12;
  if (nSides < nSidesMin) {
    nSides = nSidesMin;
    G4cout << "G4ModelingParameters::SetNoOfSides: attempt to set the"
      "\nnumber of sides per circle < " << nSidesMin
	 << "; forced to" << nSides << G4endl;
  }
  fNoOfSides = nSides;
  return fNoOfSides;
}

void G4ModelingParameters::PrintDifferences
(const G4ModelingParameters& that) const {

  if (
      (fpDefaultVisAttributes != that.fpDefaultVisAttributes) ||
      (fDrawingStyle          != that.fDrawingStyle)          ||
      (fRepStyle              != that.fRepStyle)              ||
      (fCulling               != that.fCulling)               ||
      (fCullInvisible         != that.fCullInvisible)         ||
      (fDensityCulling        != that.fDensityCulling)        ||
      (fVisibleDensity        != that.fVisibleDensity)        ||
      (fCullCovered           != that.fCullCovered)           ||
      (fNoOfSides             != that.fNoOfSides)             ||
      (fViewGeom              != that.fViewGeom)              ||
      (fViewHits              != that.fViewHits)              ||
      (fViewDigis             != that.fViewDigis))
    G4cout << "Difference in 1st batch." << G4endl;
}

std::ostream& operator << (std::ostream& os, const G4ModelingParameters& mp) {
  os << "Modeling parameters and options:";

  const G4VisAttributes* va = mp.fpDefaultVisAttributes;
  os << "\n  Default vis. attributes: ";
  if (va) os << *va;
  else os << "none";

  os << "\n  Current requested drawing style: ";
  switch (mp.fDrawingStyle) {
  case G4ModelingParameters::wf:
    os << "wireframe"; break;
  case G4ModelingParameters::hlr:
    os << "hidden line removal (hlr)"; break;
  case G4ModelingParameters::hsr:
    os << "surface (hsr)"; break;
  case G4ModelingParameters::hlhsr:
    os << "surface and edges (hlhsr)"; break;
  default: os << "unrecognised"; break;
  }

  os << "\n  Representation style for graphics reps, if needed: ";
  switch (mp.fRepStyle) {
  case G4ModelingParameters::wireframe:
    os << "wireframe"; break;
  case G4ModelingParameters::polyhedron:
    os << "polyhedron"; break;
  case G4ModelingParameters::nurbs:
    os << "nurbs"; break;
  default: os << "unrecognised"; break;
  }

  os << "\n  Culling: ";
  if (mp.fCulling) os << "on";
  else            os << "off";

  os << "\n  Culling invisible objects: ";
  if (mp.fCullInvisible) os << "on";
  else                  os << "off";

  os << "\n  Density culling: ";
  if (mp.fDensityCulling) {
    os << "on - invisible if density less than "
       << mp.fVisibleDensity / (1. * g / cm3) << " g cm^-3";
  }
  else os << "off";

  os << "\n  Culling daughters covered by opaque mothers: ";
  if (mp.fCullCovered) os << "on";
  else                os << "off";

  os << "\n  No. of sides used in circle polygon approximation: "
     << mp.fNoOfSides;

  os << "\n  View geometry: ";
  if (mp.fViewGeom) os << "true";
  else os << "false";

  os << "\n  View hits    : ";
  if (mp.fViewHits) os << "true";
  else os << "false";

  os << "\n  View digits  : ";
  if (mp.fViewDigis) os << "true";
  else os << "false";

  return os;
}

G4bool G4ModelingParameters::operator !=
(const G4ModelingParameters& mp) const {

  if (
      (*fpDefaultVisAttributes != *mp.fpDefaultVisAttributes) ||
      (fDrawingStyle           != mp.fDrawingStyle)           ||
      (fRepStyle               != mp.fRepStyle)               ||
      (fCulling                != mp.fCulling)                ||
      (fCullInvisible          != mp.fCullInvisible)          ||
      (fDensityCulling         != mp.fDensityCulling)         ||
      (fCullCovered            != mp.fCullCovered)            ||
      (fNoOfSides              != mp.fNoOfSides)              ||
      (fViewGeom               != mp.fViewGeom)               ||
      (fViewHits               != mp.fViewHits)               ||
      (fViewDigis              != mp.fViewDigis))
    return true;

  if (fDensityCulling &&
      (fVisibleDensity != mp.fVisibleDensity)) return true;

  return false;
}
