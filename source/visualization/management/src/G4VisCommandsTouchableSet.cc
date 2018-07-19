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
// $Id: G4VisCommandsTouchableSet.cc 104163 2017-05-15 06:52:42Z gcosmo $

// /vis/touchable/set commands - John Allison  8th October 2012

#include "G4VisCommandsTouchableSet.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

#include <sstream>

////////////// /vis/touchable/set/colour ///////////////////////////////////////

G4VisCommandsTouchableSet::G4VisCommandsTouchableSet()
{
  G4bool omitable;
  G4UIparameter* parameter;
  
  fpCommandSetColour = new G4UIcommand
  ("/vis/touchable/set/colour", this);
  fpCommandSetColour->SetGuidance("Set colour of current touchable.");
  fpCommandSetColour->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandSetColour->SetGuidance(ConvertToColourGuidance());
  parameter = new G4UIparameter("red", 's', omitable = true);
  parameter->SetDefaultValue("1.");
  fpCommandSetColour->SetParameter(parameter);
  parameter = new G4UIparameter("green", 'd', omitable = true);
  parameter->SetDefaultValue(1.);
  fpCommandSetColour->SetParameter(parameter);
  parameter = new G4UIparameter("blue", 'd', omitable = true);
  parameter->SetDefaultValue(1.);
  fpCommandSetColour->SetParameter(parameter);
  parameter = new G4UIparameter("opacity", 'd', omitable = true);
  parameter->SetDefaultValue(1.);
  fpCommandSetColour->SetParameter(parameter);

  fpCommandSetDaughtersInvisible = new G4UIcmdWithABool
  ("/vis/touchable/set/daughtersInvisible", this);
  fpCommandSetDaughtersInvisible->SetGuidance
  ("Daughters of current touchable invisible: true/false.");
  fpCommandSetDaughtersInvisible->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandSetDaughtersInvisible->SetParameterName("daughtersInvisible", omitable = true);
  fpCommandSetDaughtersInvisible->SetDefaultValue(true);

  fpCommandSetForceAuxEdgeVisible = new G4UIcmdWithABool
  ("/vis/touchable/set/forceAuxEdgeVisible", this);
  fpCommandSetForceAuxEdgeVisible->SetGuidance
    ("Force auxiliary (soft) edges of current touchable to be visible: true/false.");
  fpCommandSetForceAuxEdgeVisible->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandSetForceAuxEdgeVisible->SetParameterName("forceAuxEdgeVisible", omitable = true);
  fpCommandSetForceAuxEdgeVisible->SetDefaultValue(true);

  fpCommandSetLineSegmentsPerCircle = new G4UIcmdWithAnInteger
  ("/vis/touchable/set/lineSegmentsPerCircle", this);
  fpCommandSetLineSegmentsPerCircle->SetGuidance
    ("For current touchable, set number of line segments per circle, the"
     "\nprecision with which a curved line or surface is represented by a"
     "\npolygon or polyhedron, regardless of the view parameters."
     "\nNegative to pick up G4Polyhedron default value.");
  fpCommandSetLineSegmentsPerCircle->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandSetLineSegmentsPerCircle->SetParameterName("lineSegmentsPerCircle", omitable = true);
  fpCommandSetLineSegmentsPerCircle->SetDefaultValue(0);
  
  fpCommandSetForceSolid = new G4UIcmdWithABool
  ("/vis/touchable/set/forceSolid", this);
  fpCommandSetForceSolid->SetGuidance
  ("Force current touchable always to be drawn solid (surface drawing).");
  fpCommandSetForceSolid->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandSetForceSolid->SetParameterName("forceSolid", omitable = true);
  fpCommandSetForceSolid->SetDefaultValue(true);

  fpCommandSetForceWireframe = new G4UIcmdWithABool
  ("/vis/touchable/set/forceWireframe", this);
  fpCommandSetForceWireframe->SetGuidance
  ("Force current touchable always to be drawn as wireframe.");
  fpCommandSetForceWireframe->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandSetForceWireframe->SetParameterName("forceWireframe", omitable = true);
  fpCommandSetForceWireframe->SetDefaultValue(true);

  fpCommandSetLineStyle = new G4UIcmdWithAString
  ("/vis/touchable/set/lineStyle", this);
  fpCommandSetLineStyle->SetGuidance("Set line style of current touchable drawing.");
  fpCommandSetLineStyle->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandSetLineStyle->SetParameterName("lineStyle", omitable = true);
  fpCommandSetLineStyle->SetCandidates("unbroken dashed dotted");
  fpCommandSetLineStyle->SetDefaultValue("unbroken");

  fpCommandSetLineWidth = new G4UIcmdWithADouble
  ("/vis/touchable/set/lineWidth", this);
  fpCommandSetLineWidth->SetGuidance("Set line width of current touchable.");
  fpCommandSetLineWidth->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandSetLineWidth->SetParameterName("lineWidth", omitable = true);
  fpCommandSetLineWidth->SetDefaultValue(1.);

  fpCommandSetVisibility = new G4UIcmdWithABool
  ("/vis/touchable/set/visibility", this);
  fpCommandSetVisibility->SetGuidance
  ("Set visibility of current touchable: true/false.");
  fpCommandSetVisibility->SetGuidance
  ("Use \"/vis/set/touchable\" to set current touchable.");
  fpCommandSetVisibility->SetParameterName("visibility", omitable = true);
  fpCommandSetVisibility->SetDefaultValue(true);
}

G4VisCommandsTouchableSet::~G4VisCommandsTouchableSet() {
  delete fpCommandSetVisibility;
  delete fpCommandSetLineWidth;
  delete fpCommandSetLineStyle;
  delete fpCommandSetForceWireframe;
  delete fpCommandSetForceSolid;
  delete fpCommandSetLineSegmentsPerCircle;
  delete fpCommandSetForceAuxEdgeVisible;
  delete fpCommandSetDaughtersInvisible;
  delete fpCommandSetColour;
}

G4String G4VisCommandsTouchableSet::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4VisCommandsTouchableSet::SetNewValue
(G4UIcommand* command,G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: G4VisCommandsTouchableSet::SetNewValue: no current viewer."
      << G4endl;
    }
    return;
  }

  G4ViewParameters workingVP = currentViewer->GetViewParameters();
  G4VisAttributes workingVisAtts;
  
  if (command == fpCommandSetColour)
  {
    G4String redOrString;
    G4double green, blue, opacity;
    std::istringstream iss(newValue);
    iss >> redOrString >> green >> blue >> opacity;
    G4Colour colour(1,1,1,1);  // Default white and opaque.
    ConvertToColour(colour, redOrString, green, blue, opacity);    
    workingVisAtts.SetColour(colour);
    workingVP.AddVisAttributesModifier
    (G4ModelingParameters::VisAttributesModifier
     (workingVisAtts,
      G4ModelingParameters::VASColour,
      fCurrentTouchablePath));
  }
  
  else if (command == fpCommandSetDaughtersInvisible) {
    workingVisAtts.SetDaughtersInvisible(G4UIcommand::ConvertToBool(newValue));
    workingVP.AddVisAttributesModifier
    (G4ModelingParameters::VisAttributesModifier
     (workingVisAtts,
      G4ModelingParameters::VASDaughtersInvisible,
      fCurrentTouchablePath));
  }
  
  else if (command == fpCommandSetForceAuxEdgeVisible) {
    workingVisAtts.SetForceAuxEdgeVisible(G4UIcommand::ConvertToBool(newValue));
    workingVP.AddVisAttributesModifier
    (G4ModelingParameters::VisAttributesModifier
     (workingVisAtts,
      G4ModelingParameters::VASForceAuxEdgeVisible,
      fCurrentTouchablePath));
  }
  
  else if (command == fpCommandSetLineSegmentsPerCircle) {
    workingVisAtts.SetForceLineSegmentsPerCircle
    (G4UIcommand::ConvertToInt(newValue));
    workingVP.AddVisAttributesModifier
    (G4ModelingParameters::VisAttributesModifier
     (workingVisAtts,
      G4ModelingParameters::VASForceLineSegmentsPerCircle,
      fCurrentTouchablePath));
  }
  
  else if (command == fpCommandSetForceSolid) {
    workingVisAtts.SetForceSolid(G4UIcommand::ConvertToBool(newValue));
    workingVP.AddVisAttributesModifier
    (G4ModelingParameters::VisAttributesModifier
     (workingVisAtts,
      G4ModelingParameters::VASForceSolid,
      fCurrentTouchablePath));
  }
  
  else if (command == fpCommandSetForceWireframe) {
    workingVisAtts.SetForceWireframe(G4UIcommand::ConvertToBool(newValue));
    workingVP.AddVisAttributesModifier
    (G4ModelingParameters::VisAttributesModifier
     (workingVisAtts,
      G4ModelingParameters::VASForceWireframe,
      fCurrentTouchablePath));
  }
  
  else if (command == fpCommandSetLineStyle) {
    G4VisAttributes::LineStyle lineStyle = G4VisAttributes::unbroken;
    if (newValue == "dashed") {
      lineStyle = G4VisAttributes::dashed;
    } else if (newValue == "dotted") {
      lineStyle = G4VisAttributes::dotted;
    }
    // All other values are "unbroken".
    workingVisAtts.SetLineStyle(lineStyle);
    workingVP.AddVisAttributesModifier
    (G4ModelingParameters::VisAttributesModifier
     (workingVisAtts,
      G4ModelingParameters::VASLineStyle,
      fCurrentTouchablePath));
  }
  
  else if (command == fpCommandSetLineWidth) {
    workingVisAtts.SetLineWidth(G4UIcommand::ConvertToDouble(newValue));
    workingVP.AddVisAttributesModifier
    (G4ModelingParameters::VisAttributesModifier
     (workingVisAtts,
      G4ModelingParameters::VASLineWidth,
      fCurrentTouchablePath));
  }
  
  else if (command == fpCommandSetVisibility) {
    workingVisAtts.SetVisibility(G4UIcommand::ConvertToBool(newValue));
    workingVP.AddVisAttributesModifier
    (G4ModelingParameters::VisAttributesModifier
     (workingVisAtts,
      G4ModelingParameters::VASVisibility,
      fCurrentTouchablePath));
  }
  
  else {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: G4VisCommandsTouchableSet::SetNewValue: unrecognised command."
      << G4endl;
    }
    return;
  }

  SetViewParameters(currentViewer,workingVP);

  // To update all views
  G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
}
