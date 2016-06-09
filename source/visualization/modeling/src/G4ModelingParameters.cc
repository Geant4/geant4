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
// $Id$
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
  fNoOfSides             (24),
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
  const G4int  nSidesMin = 12;
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
      (fpSectionSolid          != mp.fpSectionSolid)     ||
      (fpCutawaySolid          != mp.fpCutawaySolid)     ||
      (fpEvent                 != mp.fpEvent)
      )
    return true;

  if (fDensityCulling &&
      (fVisibleDensity != mp.fVisibleDensity)) return true;

  if (VAMSNotEqual(fVisAttributesModifiers, mp.fVisAttributesModifiers))
    return true;

  return false;
}

G4bool G4ModelingParameters::PVNameCopyNoPathNotEqual
(const G4ModelingParameters::PVNameCopyNoPath& path1,
 const G4ModelingParameters::PVNameCopyNoPath& path2)
{
  if (path1.size() != path2.size()) return true;
  if (path1.empty()) return false;
  PVNameCopyNoPathConstIterator i1, i2;
  for (i1 = path1.begin(), i2 = path2.begin();
       i1 != path1.end();
       ++i1, ++i2) {
    if (i1->GetName() != i2->GetName()) return true;
    if (i1->GetCopyNo() != i2->GetCopyNo()) return true;
  }
  return false;
}

G4bool G4ModelingParameters::VAMSNotEqual
(const std::vector<G4ModelingParameters::VisAttributesModifier>& vams1,
 const std::vector<G4ModelingParameters::VisAttributesModifier>& vams2)
{
  if (vams1.size() != vams2.size()) return true;
  if (vams1.empty()) return false;
  std::vector<VisAttributesModifier>::const_iterator i1, i2;
  for (i1 = vams1.begin(), i2 = vams2.begin();
       i1 != vams1.end();
       ++i1, ++i2) {
    const PVNameCopyNoPath& vam1Path = i1->GetPVNameCopyNoPath();
    const PVNameCopyNoPath& vam2Path = i2->GetPVNameCopyNoPath();
    if (PVNameCopyNoPathNotEqual(vam1Path,vam2Path)) return true;
    if (i1->GetVisAttributesSignifier() != i2->GetVisAttributesSignifier())
      return true;
    const G4VisAttributes& vam1VisAtts = i1->GetVisAttributes();
    const G4VisAttributes& vam2VisAtts = i2->GetVisAttributes();
    const G4Colour& c1 = vam1VisAtts.GetColour();
    const G4Colour& c2 = vam2VisAtts.GetColour();
    switch (i1->GetVisAttributesSignifier()) {
      case G4ModelingParameters::VASVisibility:
        if (vam1VisAtts.IsVisible() != vam2VisAtts.IsVisible())
          return true;
        break;
      case G4ModelingParameters::VASDaughtersInvisible:
        if (vam1VisAtts.IsDaughtersInvisible() !=
            vam2VisAtts.IsDaughtersInvisible())
          return true;
        break;
      case G4ModelingParameters::VASColour:
        if (c1 != c2)
          return true;
        break;
      case G4ModelingParameters::VASLineStyle:
        if (vam1VisAtts.GetLineStyle() != vam2VisAtts.GetLineStyle())
          return true;
        break;
      case G4ModelingParameters::VASLineWidth:
        if (vam1VisAtts.GetLineWidth() != vam2VisAtts.GetLineWidth())
          return true;
        break;
      case G4ModelingParameters::VASForceWireframe:
        if (vam1VisAtts.GetForcedDrawingStyle() !=
            vam2VisAtts.GetForcedDrawingStyle())
          return true;
        break;
      case G4ModelingParameters::VASForceSolid:
        if (vam1VisAtts.GetForcedDrawingStyle() !=
            vam2VisAtts.GetForcedDrawingStyle())
          return true;
        break;
      case G4ModelingParameters::VASForceAuxEdgeVisible:
        if (vam1VisAtts.IsForceAuxEdgeVisible() !=
            vam2VisAtts.IsForceAuxEdgeVisible())
          return true;
        break;
      case G4ModelingParameters::VASForceLineSegmentsPerCircle:
        if (vam1VisAtts.GetForcedLineSegmentsPerCircle() !=
            vam2VisAtts.GetForcedLineSegmentsPerCircle())
          return true;
        break;
    }
  }
  return false;
}

std::ostream& operator <<
(std::ostream& os, const G4ModelingParameters::PVNameCopyNoPath& path)
{
//  os << "Touchable path: physical-volume-name:copy-number pairs:\n  ";
  G4ModelingParameters::PVNameCopyNoPathConstIterator i;
  for (i = path.begin(); i != path.end(); ++i) {
    if (i != path.begin()) {
      os << ", ";
    }
    os << i->GetName() << ':' << i->GetCopyNo();
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
        os << " forceAuxEdgeVisible ";
        if (vamVisAtts.IsForceAuxEdgeVisible()) {
          os << "true";
        } else {
          os << "false";
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

