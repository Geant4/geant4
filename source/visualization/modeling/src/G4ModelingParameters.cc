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
// $Id: G4ModelingParameters.cc 103627 2017-04-19 13:30:26Z gcosmo $
//
// 
// John Allison  31st December 1997.
// Parameters associated with the modeling of GEANT4 objects.

#include "G4ModelingParameters.hh"

#include "G4ios.hh"
#include "G4VisAttributes.hh"
#include "G4ExceptionSeverity.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"

G4ModelingParameters::G4ModelingParameters ():
  fWarning               (true),
  fpDefaultVisAttributes (0),
  fDrawingStyle          (wf),
  fCulling               (false),
  fCullInvisible         (false),
  fDensityCulling        (false),
  fVisibleDensity        (0.01 * g / cm3),
  fCullCovered           (false),
  fExplodeFactor         (1.),
  fNoOfSides             (72),
  fpSectionSolid         (0),
  fpCutawaySolid         (0),
  fpEvent                (0)
{}

G4ModelingParameters::G4ModelingParameters
(const G4VisAttributes* pDefaultVisAttributes,
 G4ModelingParameters::DrawingStyle drawingStyle,
 G4bool isCulling,
 G4bool isCullingInvisible,
 G4bool isDensityCulling,
 G4double visibleDensity,
 G4bool isCullingCovered,
 G4int noOfSides
 ):
  fWarning        (true),
  fpDefaultVisAttributes (pDefaultVisAttributes),
  fDrawingStyle   (drawingStyle),
  fCulling        (isCulling),
  fCullInvisible  (isCullingInvisible),
  fDensityCulling (isDensityCulling),
  fVisibleDensity (visibleDensity),
  fCullCovered    (isCullingCovered),
  fExplodeFactor  (1.),
  fNoOfSides      (noOfSides),
  fpSectionSolid  (0),
  fpCutawaySolid  (0),
  fpEvent         (0)
{}

G4ModelingParameters::~G4ModelingParameters ()
{
  delete fpSectionSolid;
  delete fpCutawaySolid;
}

G4ModelingParameters::VisAttributesModifier::VisAttributesModifier
(const G4VisAttributes& visAtts,
 G4ModelingParameters::VisAttributesSignifier signifier,
 const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& path):
fVisAtts(visAtts), fSignifier(signifier)
{
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
  typedef PVPath::const_iterator PVPathConstIterator;
  PVPathConstIterator i;
  for (i = path.begin();
       i != path.end();
       ++i) {
    fPVNameCopyNoPath.push_back
    (PVNameCopyNo
     (i->GetPhysicalVolume()->GetName(),
      i->GetCopyNo()));
  }
}

void G4ModelingParameters::SetVisibleDensity (G4double visibleDensity) {
  const G4double reasonableMaximum = 10.0 * g / cm3;
  if (visibleDensity < 0 && fWarning) {
    G4cout << "G4ModelingParameters::SetVisibleDensity: attempt to set negative "
      "density - ignored." << G4endl;
  }
  else {
    if (fVisibleDensity > reasonableMaximum && fWarning) {
      G4cout << "G4ModelingParameters::SetVisibleDensity: density > "
	   << reasonableMaximum
	   << " g / cm3 - did you mean this?"
	   << G4endl;
    }
    fVisibleDensity = visibleDensity;
  }
}

G4int G4ModelingParameters::SetNoOfSides (G4int nSides) {
  const G4int  nSidesMin = fpDefaultVisAttributes->GetMinLineSegmentsPerCircle();
  if (nSides < nSidesMin) {
    nSides = nSidesMin;
    if (fWarning)
      G4cout << "G4ModelingParameters::SetNoOfSides: attempt to set the"
	"\nnumber of sides per circle < " << nSidesMin
	     << "; forced to" << nSides << G4endl;
  }
  fNoOfSides = nSides;
  return fNoOfSides;
}

void G4ModelingParameters::SetSectionSolid
(G4VSolid* pSectionSolid) {
  delete fpSectionSolid;
  fpSectionSolid = pSectionSolid;
}

void G4ModelingParameters::SetCutawaySolid
(G4VSolid* pCutawaySolid) {
  delete fpCutawaySolid;
  fpCutawaySolid = pCutawaySolid;
}

std::ostream& operator << (std::ostream& os, const G4ModelingParameters& mp)
{
  os << "Modeling parameters (warning ";
  if (mp.fWarning) os << "true";
  else os << "false";
  os << "):";

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

  os << "\n  Explode factor: " << mp.fExplodeFactor
     << " about centre: " << mp.fExplodeCentre;

  os << "\n  No. of sides used in circle polygon approximation: "
     << mp.fNoOfSides;

  os << "\n  Section (DCUT) shape (G4VSolid) pointer: ";
  if (!mp.fpSectionSolid) os << "non-";
  os << "null";

  os << "\n  Cutaway (DCUT) shape (G4VSolid) pointer: ";
  if (!mp.fpCutawaySolid) os << "non-";
  os << "null";

  os << "\n  Event pointer: " << mp.fpEvent;

  os << "\n  Vis attributes modifiers: ";
  const std::vector<G4ModelingParameters::VisAttributesModifier>& vams =
  mp.fVisAttributesModifiers;
  if (vams.empty()) {
    os << "None";
  } else {
    os << vams;
  }

  return os;
}

G4bool G4ModelingParameters::operator !=
(const G4ModelingParameters& mp) const {

  if (
      (fWarning                != mp.fWarning)                ||
      (*fpDefaultVisAttributes != *mp.fpDefaultVisAttributes) ||
      (fCulling                != mp.fCulling)                ||
      (fCullInvisible          != mp.fCullInvisible)          ||
      (fDensityCulling         != mp.fDensityCulling)         ||
      (fCullCovered            != mp.fCullCovered)            ||
      (fExplodeFactor          != mp.fExplodeFactor)          ||
      (fExplodeCentre          != mp.fExplodeCentre)          ||
      (fNoOfSides              != mp.fNoOfSides)              ||
      (fpSectionSolid          != mp.fpSectionSolid)          ||
      (fpCutawaySolid          != mp.fpCutawaySolid)          ||
      (fpEvent                 != mp.fpEvent)
      )
    return true;

  if (fDensityCulling &&
      (fVisibleDensity != mp.fVisibleDensity)) return true;

  if (fVisAttributesModifiers != mp.fVisAttributesModifiers)
    return true;

  return false;
}

G4bool G4ModelingParameters::VisAttributesModifier::operator!=
(const G4ModelingParameters::VisAttributesModifier& rhs) const
{
  if (fSignifier != rhs.fSignifier) return true;
  if (fPVNameCopyNoPath != rhs.fPVNameCopyNoPath) return true;
  switch (fSignifier) {
    case G4ModelingParameters::VASVisibility:
      if (fVisAtts.IsVisible() != rhs.fVisAtts.IsVisible())
        return true;
      break;
    case G4ModelingParameters::VASDaughtersInvisible:
      if (fVisAtts.IsDaughtersInvisible() !=
          rhs.fVisAtts.IsDaughtersInvisible())
        return true;
      break;
    case G4ModelingParameters::VASColour:
      if (fVisAtts.GetColour() != rhs.fVisAtts.GetColour())
        return true;
      break;
    case G4ModelingParameters::VASLineStyle:
      if (fVisAtts.GetLineStyle() != rhs.fVisAtts.GetLineStyle())
        return true;
      break;
    case G4ModelingParameters::VASLineWidth:
      if (fVisAtts.GetLineWidth() != rhs.fVisAtts.GetLineWidth())
        return true;
      break;
    case G4ModelingParameters::VASForceWireframe:
      if (fVisAtts.GetForcedDrawingStyle() !=
          rhs.fVisAtts.GetForcedDrawingStyle())
        return true;
      break;
    case G4ModelingParameters::VASForceSolid:
      if (fVisAtts.GetForcedDrawingStyle() !=
          rhs.fVisAtts.GetForcedDrawingStyle())
        return true;
      break;
    case G4ModelingParameters::VASForceAuxEdgeVisible:
      if (fVisAtts.IsForceAuxEdgeVisible() !=
          rhs.fVisAtts.IsForceAuxEdgeVisible() ||
          fVisAtts.IsForcedAuxEdgeVisible() !=
          rhs.fVisAtts.IsForcedAuxEdgeVisible())
        return true;
      break;
    case G4ModelingParameters::VASForceLineSegmentsPerCircle:
      if (fVisAtts.GetForcedLineSegmentsPerCircle() !=
          rhs.fVisAtts.GetForcedLineSegmentsPerCircle())
        return true;
      break;
  }
  return false;
}

G4bool G4ModelingParameters::PVNameCopyNo::operator!=
(const G4ModelingParameters::PVNameCopyNo& rhs) const
{
  if (fName != rhs.fName) return true;
  if (fCopyNo != rhs.fCopyNo) return true;
  return false;
}

std::ostream& operator <<
(std::ostream& os, const G4ModelingParameters::PVNameCopyNoPath& path)
{
  os << "Touchable path: physical-volume-name:copy-number pairs:\n  ";
  G4ModelingParameters::PVNameCopyNoPathConstIterator i;
  for (i = path.begin(); i != path.end(); ++i) {
    if (i != path.begin()) {
      os << ',';
    }
    os << i->GetName() << ':' << i->GetCopyNo();
  }
  return os;
}

const G4String& G4ModelingParameters::PVPointerCopyNo::GetName() const
{
  return fpPV->GetName();
}

G4bool G4ModelingParameters::PVPointerCopyNo::operator!=
(const G4ModelingParameters::PVPointerCopyNo& rhs) const
{
  if (fpPV != rhs.fpPV) return true;
  if (fCopyNo != rhs.fCopyNo) return true;
  return false;
}

std::ostream& operator <<
(std::ostream& os, const G4ModelingParameters::PVPointerCopyNoPath& path)
{
  os << "Touchable path: physical-volume-pointer:copy-number pairs:\n  ";
  G4ModelingParameters::PVPointerCopyNoPathConstIterator i;
  for (i = path.begin(); i != path.end(); ++i) {
    if (i != path.begin()) {
      os << ',';
    }
    os << '(' << (void*)(i->GetPVPointer()) << ')' << i->GetName() << ':' << i->GetCopyNo();
  }
  return os;
}

std::ostream& operator <<
(std::ostream& os,
 const std::vector<G4ModelingParameters::VisAttributesModifier>& vams)
{
  std::vector<G4ModelingParameters::VisAttributesModifier>::const_iterator
  iModifier;
  for (iModifier = vams.begin();
       iModifier != vams.end();
       ++iModifier) {
    const G4ModelingParameters::PVNameCopyNoPath& vamPath =
    iModifier->GetPVNameCopyNoPath();
    os << '\n' << vamPath;
    const G4VisAttributes& vamVisAtts = iModifier->GetVisAttributes();
    const G4Colour& c = vamVisAtts.GetColour();
    switch (iModifier->GetVisAttributesSignifier()) {
      case G4ModelingParameters::VASVisibility:
        os << " visibility ";
        if (vamVisAtts.IsVisible()) {
          os << "true";
        } else {
          os << "false";
        }
        break;
      case G4ModelingParameters::VASDaughtersInvisible:
        os << " daughtersInvisible ";
        if (vamVisAtts.IsDaughtersInvisible()) {
          os << "true";
        } else {
          os << "false";
        }
        break;
      case G4ModelingParameters::VASColour:
        os << " colour " << c;
        break;
      case G4ModelingParameters::VASLineStyle:
        os << " lineStyle ";
        switch (vamVisAtts.GetLineStyle()) {
          case G4VisAttributes::unbroken:
            os << "unbroken";
            break;
          case G4VisAttributes::dashed:
            os << "dashed";
            break;
          case G4VisAttributes::dotted:
            os << "dotted";
        }
        break;
      case G4ModelingParameters::VASLineWidth:
        os << " lineWidth "
        << vamVisAtts.GetLineWidth();
        break;
      case G4ModelingParameters::VASForceWireframe:
        if (vamVisAtts.GetForcedDrawingStyle() == G4VisAttributes::wireframe) {
          os << " forceWireframe ";
          if (vamVisAtts.IsForceDrawingStyle()) {
            os << "true";
          } else {
            os << "false";
          }
        }
        break;
      case G4ModelingParameters::VASForceSolid:
        if (vamVisAtts.GetForcedDrawingStyle() == G4VisAttributes::solid) {
          os << " forceSolid ";
          if (vamVisAtts.IsForceDrawingStyle()) {
            os << "true";
          } else {
            os << "false";
          }
        }
        break;
      case G4ModelingParameters::VASForceAuxEdgeVisible:
        os << " forceAuxEdgeVisible: ";
        if (!vamVisAtts.IsForceDrawingStyle()) {
          os << "not ";
        }
        os << " forced";
        if (vamVisAtts.IsForceAuxEdgeVisible()) {
          os << ": ";
          if (vamVisAtts.IsForcedAuxEdgeVisible()) {
            os << "true";
          } else {
            os << "false";
          }
        }
        break;
      case G4ModelingParameters::VASForceLineSegmentsPerCircle:
        os << " lineSegmentsPerCircle "
        << vamVisAtts.GetForcedLineSegmentsPerCircle();
        break;
    }
  }

  return os;
}

