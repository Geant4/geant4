// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ModelingParameters.cc,v 1.2 1999-01-10 13:25:50 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  31st December 1997.
// Parameters associated with the modeling of GEANT4 objects.

#include "G4ModelingParameters.hh"

#include "G4ios.hh"
#include "G4VisAttributes.hh"

G4ModelingParameters::G4ModelingParameters ():
  fpDefaultVisAttributes (0),
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
      "density - ignored." << endl;
  }
  else {
    if (fVisibleDensity > reasonableMaximum) {
      G4cout << "G4ModelingParameters::SetVisibleDensity: density > "
	   << reasonableMaximum
	   << " g / cm3 - did you mean this?"
	   << endl;
    }
    fVisibleDensity = visibleDensity;
  }
}

void G4ModelingParameters::SetNoOfSides (G4int nSides) {
  const G4int  nSidesMin = 3;
  if (nSides < nSidesMin) {
    nSides = nSidesMin;
    G4cout << "G4ModelingParameters::SetNoOfSides: attempt to set the"
      "\nnumber of sides per circle < " << nSidesMin
	 << "; forced to" << nSides << endl;
  }
  fNoOfSides = nSides;
}

void G4ModelingParameters::PrintDifferences
(const G4ModelingParameters& that) const {

  if (
      (fpDefaultVisAttributes != that.fpDefaultVisAttributes) ||
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
    G4cout << "Difference in 1st batch." << endl;
}

ostream& operator << (ostream& os, const G4ModelingParameters& mp) {
  os << "Modeling parameters and options:";

  os << "\n  Default vis. attributes: " << *mp.fpDefaultVisAttributes;

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

G4bool operator != (const G4ModelingParameters& mp1,
		    const G4ModelingParameters& mp2) {

  if (
      (*mp1.fpDefaultVisAttributes != *mp2.fpDefaultVisAttributes) ||
      (mp1.fRepStyle               != mp2.fRepStyle)               ||
      (mp1.fCulling                != mp2.fCulling)                ||
      (mp1.fCullInvisible          != mp2.fCullInvisible)          ||
      (mp1.fDensityCulling         != mp2.fDensityCulling)         ||
      (mp1.fCullCovered            != mp2.fCullCovered)            ||
      (mp1.fNoOfSides              != mp2.fNoOfSides)              ||
      (mp1.fViewGeom               != mp2.fViewGeom)               ||
      (mp1.fViewHits               != mp2.fViewHits)               ||
      (mp1.fViewDigis              != mp2.fViewDigis))
    return true;

  if (mp1.fDensityCulling &&
      (mp1.fVisibleDensity != mp2.fVisibleDensity)) return true;

  return false;
}
