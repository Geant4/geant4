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
// $Id: G4VisCommandsGeometrySet.cc 104123 2017-05-11 13:51:32Z gcosmo $

// /vis/geometry commands - John Allison  31st January 2006

#include "G4VisCommandsGeometrySet.hh"

#include "G4UIcommand.hh"
#include "G4VisManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4UImanager.hh"

#include <sstream>

void G4VVisCommandGeometrySet::Set
(G4String requestedName,
 const G4VVisCommandGeometrySetFunction& setFunction,
 G4int requestedDepth)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  G4bool found = false;
  for (size_t iLV = 0; iLV < pLVStore->size(); iLV++ ) {
    G4LogicalVolume* pLV = (*pLVStore)[iLV];
    const G4String& logVolName = pLV->GetName();
    if (logVolName == requestedName) found = true;
    if (requestedName == "all" || logVolName == requestedName) {
      SetLVVisAtts(pLV, setFunction, 0, requestedDepth);
    }
  }
  if (requestedName != "all" && !found) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Logical volume \"" << requestedName
	     << "\" not found in logical volume store." << G4endl;
    }
    return;
  }
  if (fpVisManager->GetCurrentViewer()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
  }
}

void G4VVisCommandGeometrySet::SetLVVisAtts
(G4LogicalVolume* pLV,
 const G4VVisCommandGeometrySetFunction& setFunction,
 G4int depth, G4int requestedDepth)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  const G4VisAttributes* oldVisAtts = pLV->GetVisAttributes();
  fVisAttsMap.insert(std::make_pair(pLV,oldVisAtts));  // Store old vis atts.
  G4VisAttributes* newVisAtts = new G4VisAttributes;   // Memory leak!
  if (oldVisAtts) {
    *newVisAtts = *oldVisAtts;
  }
  setFunction(newVisAtts);  // Sets whatever attribute determined by
			    // function object.
  pLV->SetVisAttributes(newVisAtts);
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "\nLogical Volume \"" << pLV->GetName()
	   << "\": setting vis attributes:";
    if (oldVisAtts) {
      G4cout << "\nwas: " << *oldVisAtts;
    } else {
      G4cout << "\n(no old attributes)";
    }
    G4cout << "\nnow: " << *newVisAtts
	   << G4endl;
  }
  if (requestedDepth < 0 || depth < requestedDepth) {
    G4int nDaughters = pLV->GetNoDaughters();
    for (G4int i = 0; i < nDaughters; ++i) {
      SetLVVisAtts(pLV->GetDaughter(i)->GetLogicalVolume(),
		   setFunction, ++depth, requestedDepth);
    }
  }
}

////////////// /vis/geometry/set/colour ///////////////////////////////////////

G4VisCommandGeometrySetColour::G4VisCommandGeometrySetColour()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/geometry/set/colour", this);
  fpCommand->SetGuidance("Sets colour of logical volume(s).");
  fpCommand->SetGuidance("\"all\" sets all logical volumes.");
  fpCommand->SetGuidance
    ("Optionally propagates down hierarchy to given depth.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = true);
  parameter->SetDefaultValue("all");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth", 'd', omitable = true);
  parameter->SetDefaultValue(0);
  parameter->SetGuidance
    ("Depth of propagation (-1 means unlimited depth).");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("red", 's', omitable = true);
  parameter->SetDefaultValue("1.");
  parameter->SetGuidance
    ("Red component or a string, e.g., \"blue\", in which case succeeding colour components are ignored.");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("green", 'd', omitable = true);
  parameter->SetDefaultValue(1.);
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("blue", 'd', omitable = true);
  parameter->SetDefaultValue(1.);
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("opacity", 'd', omitable = true);
  parameter->SetDefaultValue(1.);
  fpCommand->SetParameter(parameter);
}

G4VisCommandGeometrySetColour::~G4VisCommandGeometrySetColour()
{
  delete fpCommand;
}

G4String G4VisCommandGeometrySetColour::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometrySetColour::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4String name, redOrString;
  G4int requestedDepth;
  G4double green, blue, opacity;
  std::istringstream iss(newValue);
  iss >> name >> requestedDepth >> redOrString >> green >> blue >> opacity;
  G4Colour colour(1,1,1,1);  // Default white and opaque.
  ConvertToColour(colour, redOrString, green, blue, opacity);
  G4VisCommandGeometrySetColourFunction setColour(colour);
  Set(name, setColour, requestedDepth);
}

////////////// /vis/geometry/set/daughtersInvisible //////////////////////

G4VisCommandGeometrySetDaughtersInvisible::G4VisCommandGeometrySetDaughtersInvisible()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/geometry/set/daughtersInvisible", this);
  fpCommand->SetGuidance("Makes daughters of logical volume(s) invisible.");
  fpCommand->SetGuidance("\"all\" sets all logical volumes.");
  fpCommand->SetGuidance
    ("Optionally propagates down hierarchy to given depth.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = true);
  parameter->SetDefaultValue("all");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth", 'd', omitable = true);
  parameter->SetDefaultValue(0);
  parameter->SetGuidance
    ("Depth of propagation (-1 means unlimited depth).");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("daughtersInvisible", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  fpCommand->SetParameter(parameter);
}

G4VisCommandGeometrySetDaughtersInvisible::~G4VisCommandGeometrySetDaughtersInvisible()
{
  delete fpCommand;
}

G4String
G4VisCommandGeometrySetDaughtersInvisible::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometrySetDaughtersInvisible::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4String name;
  G4int requestedDepth;
  G4String daughtersInvisibleString;
  std::istringstream iss(newValue);
  iss >> name >> requestedDepth >> daughtersInvisibleString;
  G4bool daughtersInvisible =
    G4UIcommand::ConvertToBool(daughtersInvisibleString);

  if (requestedDepth !=0) {
    requestedDepth = 0;
    if (fpVisManager->GetVerbosity() >= G4VisManager::warnings) {
      G4cout << "Recursive application suppressed for this attribute."
	     << G4endl;
    }
  }

  G4VisCommandGeometrySetDaughtersInvisibleFunction
    setDaughtersInvisible(daughtersInvisible);
  Set(name, setDaughtersInvisible, requestedDepth);

  G4VViewer* pViewer = fpVisManager->GetCurrentViewer();
  if (pViewer) {
    const G4ViewParameters& viewParams = pViewer->GetViewParameters();
    if (fpVisManager->GetVerbosity() >= G4VisManager::warnings) {
      if (!viewParams.IsCulling()) {
	G4cout <<
	  "Culling must be on - \"/vis/viewer/set/culling global true\" - to see effect."
	       << G4endl;
      }
    }
  }
}

////////////// /vis/geometry/set/forceAuxEdgeVisible /////////////////////////

G4VisCommandGeometrySetForceAuxEdgeVisible::G4VisCommandGeometrySetForceAuxEdgeVisible()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/geometry/set/forceAuxEdgeVisible", this);
  fpCommand->SetGuidance
    ("Forces auxiliary (soft) edges of logical volume(s) to be visible,"
    "\nregardless of the view parameters.");
  fpCommand->SetGuidance("\"all\" sets all logical volumes.");
  fpCommand->SetGuidance
    ("Optionally propagates down hierarchy to given depth.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = true);
  parameter->SetDefaultValue("all");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth", 'd', omitable = true);
  parameter->SetDefaultValue(0);
  parameter->SetGuidance
    ("Depth of propagation (-1 means unlimited depth).");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("forceAuxEdgeVisible", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  fpCommand->SetParameter(parameter);
}

G4VisCommandGeometrySetForceAuxEdgeVisible::~G4VisCommandGeometrySetForceAuxEdgeVisible()
{
  delete fpCommand;
}

G4String
G4VisCommandGeometrySetForceAuxEdgeVisible::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometrySetForceAuxEdgeVisible::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4String name;
  G4int requestedDepth;
  G4String forceAuxEdgeVisibleString;
  std::istringstream iss(newValue);
  iss >> name >> requestedDepth >> forceAuxEdgeVisibleString;
  G4bool forceAuxEdgeVisible =
    G4UIcommand::ConvertToBool(forceAuxEdgeVisibleString);;

  G4VisCommandGeometrySetForceAuxEdgeVisibleFunction
    setForceAuxEdgeVisible(forceAuxEdgeVisible);
  Set(name, setForceAuxEdgeVisible, requestedDepth);
}

////////////// /vis/geometry/set/forceLineSegmentsPerCircle /////////////////////////

G4VisCommandGeometrySetForceLineSegmentsPerCircle::G4VisCommandGeometrySetForceLineSegmentsPerCircle()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/geometry/set/forceLineSegmentsPerCircle", this);
  fpCommand->SetGuidance
    ("Forces number of line segments per circle, the precision with which a"
     "\ncurved line or surface is represented by a polygon or polyhedron,"
     "\nregardless of the view parameters.");
  fpCommand->SetGuidance("\"all\" sets all logical volumes.");
  fpCommand->SetGuidance
    ("Optionally propagates down hierarchy to given depth.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = true);
  parameter->SetDefaultValue("all");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth", 'd', omitable = true);
  parameter->SetDefaultValue(0);
  parameter->SetGuidance
    ("Depth of propagation (-1 means unlimited depth).");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("lineSegmentsPerCircle", 'd', omitable = true);
  parameter->SetGuidance
    ("<= 0 means not forced, i.e., under control of viewer.");
  parameter->SetDefaultValue(0);
  fpCommand->SetParameter(parameter);
}

G4VisCommandGeometrySetForceLineSegmentsPerCircle::~G4VisCommandGeometrySetForceLineSegmentsPerCircle()
{
  delete fpCommand;
}

G4String
G4VisCommandGeometrySetForceLineSegmentsPerCircle::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometrySetForceLineSegmentsPerCircle::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4String name;
  G4int requestedDepth;
  G4int lineSegmentsPerCircle;
  std::istringstream iss(newValue);
  iss >> name >> requestedDepth >> lineSegmentsPerCircle;

  G4VisCommandGeometrySetForceLineSegmentsPerCircleFunction setForceLineSegmentsPerCircle(lineSegmentsPerCircle);
  Set(name, setForceLineSegmentsPerCircle, requestedDepth);
}

////////////// /vis/geometry/set/forceSolid /////////////////////////

G4VisCommandGeometrySetForceSolid::G4VisCommandGeometrySetForceSolid()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/geometry/set/forceSolid", this);
  fpCommand->SetGuidance
   ("Forces logical volume(s) always to be drawn solid (surface drawing),"
    "\nregardless of the view parameters.");
  fpCommand->SetGuidance("\"all\" sets all logical volumes.");
  fpCommand->SetGuidance
    ("Optionally propagates down hierarchy to given depth.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = true);
  parameter->SetDefaultValue("all");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth", 'd', omitable = true);
  parameter->SetDefaultValue(0);
  parameter->SetGuidance
    ("Depth of propagation (-1 means unlimited depth).");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("forceSolid", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  fpCommand->SetParameter(parameter);
}

G4VisCommandGeometrySetForceSolid::~G4VisCommandGeometrySetForceSolid()
{
  delete fpCommand;
}

G4String
G4VisCommandGeometrySetForceSolid::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometrySetForceSolid::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4String name;
  G4int requestedDepth;
  G4String forceSolidString;
  std::istringstream iss(newValue);
  iss >> name >> requestedDepth >> forceSolidString;
  G4bool forceSolid = G4UIcommand::ConvertToBool(forceSolidString);

  G4VisCommandGeometrySetForceSolidFunction setForceSolid(forceSolid);
  Set(name, setForceSolid, requestedDepth);
}

////////////// /vis/geometry/set/forceWireframe /////////////////////////

G4VisCommandGeometrySetForceWireframe::G4VisCommandGeometrySetForceWireframe()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/geometry/set/forceWireframe", this);
  fpCommand->SetGuidance
   ("Forces logical volume(s) always to be drawn as wireframe,"
    "\nregardless of the view parameters.");
  fpCommand->SetGuidance("\"all\" sets all logical volumes.");
  fpCommand->SetGuidance
    ("Optionally propagates down hierarchy to given depth.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = true);
  parameter->SetDefaultValue("all");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth", 'd', omitable = true);
  parameter->SetDefaultValue(0);
  parameter->SetGuidance
    ("Depth of propagation (-1 means unlimited depth).");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("forceWireframe", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  fpCommand->SetParameter(parameter);
}

G4VisCommandGeometrySetForceWireframe::~G4VisCommandGeometrySetForceWireframe()
{
  delete fpCommand;
}

G4String
G4VisCommandGeometrySetForceWireframe::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometrySetForceWireframe::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4String name;
  G4int requestedDepth;
  G4String forceWireframeString;
  std::istringstream iss(newValue);
  iss >> name >> requestedDepth >> forceWireframeString;
  G4bool forceWireframe = G4UIcommand::ConvertToBool(forceWireframeString);

  G4VisCommandGeometrySetForceWireframeFunction
    setForceWireframe(forceWireframe);
  Set(name, setForceWireframe, requestedDepth);
}

////////////// /vis/geometry/set/lineStyle /////////////////////////////////

G4VisCommandGeometrySetLineStyle::G4VisCommandGeometrySetLineStyle()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/geometry/set/lineStyle", this);
  fpCommand->SetGuidance("Sets line style of logical volume(s) drawing.");
  fpCommand->SetGuidance("\"all\" sets all logical volumes.");
  fpCommand->SetGuidance
    ("Optionally propagates down hierarchy to given depth.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = true);
  parameter->SetDefaultValue("all");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth", 'd', omitable = true);
  parameter->SetDefaultValue(0);
  parameter->SetGuidance
    ("Depth of propagation (-1 means unlimited depth).");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("lineStyle", 's', omitable = true);
  parameter->SetParameterCandidates("unbroken dashed dotted");
  parameter->SetDefaultValue("unbroken");
  fpCommand->SetParameter(parameter);
}

G4VisCommandGeometrySetLineStyle::~G4VisCommandGeometrySetLineStyle()
{
  delete fpCommand;
}

G4String
G4VisCommandGeometrySetLineStyle::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometrySetLineStyle::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4String name, lineStyleString;
  G4int requestedDepth;
  std::istringstream iss(newValue);
  iss >> name >> requestedDepth >> lineStyleString;
  G4VisAttributes::LineStyle lineStyle = G4VisAttributes::unbroken;
  if (lineStyleString == "unbroken") lineStyle = G4VisAttributes::unbroken;
  if (lineStyleString == "dashed") lineStyle = G4VisAttributes::dashed;
  if (lineStyleString == "dotted") lineStyle = G4VisAttributes::dotted;

  G4VisCommandGeometrySetLineStyleFunction setLineStyle(lineStyle);
  Set(name, setLineStyle, requestedDepth);
}

////////////// /vis/geometry/set/lineWidth /////////////////////////////////

G4VisCommandGeometrySetLineWidth::G4VisCommandGeometrySetLineWidth()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/geometry/set/lineWidth", this);
  fpCommand->SetGuidance("Sets line width of logical volume(s) drawing.");
  fpCommand->SetGuidance("\"all\" sets all logical volumes.");
  fpCommand->SetGuidance
    ("Optionally propagates down hierarchy to given depth.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = true);
  parameter->SetDefaultValue("all");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth", 'd', omitable = true);
  parameter->SetDefaultValue(0);
  parameter->SetGuidance
    ("Depth of propagation (-1 means unlimited depth).");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("lineWidth", 'd', omitable = true);
  parameter->SetDefaultValue(1.);
  fpCommand->SetParameter(parameter);
}

G4VisCommandGeometrySetLineWidth::~G4VisCommandGeometrySetLineWidth()
{
  delete fpCommand;
}

G4String
G4VisCommandGeometrySetLineWidth::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometrySetLineWidth::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4String name;
  G4int requestedDepth;
  G4double lineWidth;
  std::istringstream iss(newValue);
  iss >> name >> requestedDepth >> lineWidth;

  G4VisCommandGeometrySetLineWidthFunction setLineWidth(lineWidth);
  Set(name, setLineWidth, requestedDepth);
}

////////////// /vis/geometry/set/visibility ///////////////////////////////////////

G4VisCommandGeometrySetVisibility::G4VisCommandGeometrySetVisibility()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/geometry/set/visibility", this);
  fpCommand->SetGuidance("Sets visibility of logical volume(s).");
  fpCommand->SetGuidance("\"all\" sets all logical volumes.");
  fpCommand->SetGuidance
    ("Optionally propagates down hierarchy to given depth.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = true);
  parameter->SetDefaultValue("all");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth", 'd', omitable = true);
  parameter->SetDefaultValue(0);
  parameter->SetGuidance
    ("Depth of propagation (-1 means unlimited depth).");
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("visibility", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  fpCommand->SetParameter(parameter);
}

G4VisCommandGeometrySetVisibility::~G4VisCommandGeometrySetVisibility()
{
  delete fpCommand;
}

G4String G4VisCommandGeometrySetVisibility::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandGeometrySetVisibility::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4String name;
  G4int requestedDepth;
  G4String visibilityString;
  std::istringstream iss(newValue);
  iss >> name >> requestedDepth >> visibilityString;
  G4bool visibility = G4UIcommand::ConvertToBool(visibilityString);

  G4VisCommandGeometrySetVisibilityFunction setVisibility(visibility);
  Set(name, setVisibility, requestedDepth);

  G4VViewer* pViewer = fpVisManager->GetCurrentViewer();
  if (pViewer) {
    const G4ViewParameters& viewParams = pViewer->GetViewParameters();
    if (fpVisManager->GetVerbosity() >= G4VisManager::warnings) {
      if (!viewParams.IsCulling() ||
	  !viewParams.IsCullingInvisible()) {
	G4cout <<
	  "Culling must be on - \"/vis/viewer/set/culling global true\" and"
	  "\n  \"/vis/viewer/set/culling invisible true\" - to see effect."
	       << G4endl;
      }
    }
  }
}

void G4VisCommandGeometrySetVisibility::SetNewValueOnLV
(G4LogicalVolume* pLV, G4int requestedDepth,G4bool visibility)
{
  if (!pLV) return;
  G4VisCommandGeometrySetVisibilityFunction setVisibility(visibility);
  SetLVVisAtts(pLV, setVisibility, 0, requestedDepth);

  G4VViewer* pViewer = fpVisManager->GetCurrentViewer();
  if (pViewer) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
    const G4ViewParameters& viewParams = pViewer->GetViewParameters();
    if (fpVisManager->GetVerbosity() >= G4VisManager::warnings) {
      if (!viewParams.IsCulling() ||
	  !viewParams.IsCullingInvisible()) {
	G4cout <<
	  "Culling must be on - \"/vis/viewer/set/culling global true\" and"
	  "\n  \"/vis/viewer/set/culling invisible true\" - to see effect."
	       << G4endl;
      }
    }
  }
}
