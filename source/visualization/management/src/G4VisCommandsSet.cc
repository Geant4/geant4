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

// /vis/set - John Allison  21st March 2012
// Set quantities for use in appropriate future commands.

#include "G4VisCommandsSet.hh"

#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4TouchablePropertiesScene.hh"
#include "G4TransportationManager.hh"
#include "G4BoundingExtentScene.hh"
#include <sstream>

////////////// /vis/set/arrow3DLineSegmentsPerCircle ////////////////////////////////////

G4VisCommandSetArrow3DLineSegmentsPerCircle::G4VisCommandSetArrow3DLineSegmentsPerCircle ()
{
  G4bool omitable;
  fpCommand = new G4UIcmdWithAnInteger("/vis/set/arrow3DLineSegmentsPerCircle", this);
  fpCommand->SetGuidance
  ("Defines number of line segments per circle for drawing 3D arrows"
   " for future \"/vis/scene/add/\" commands.");
  fpCommand->SetParameterName ("number", omitable = true);
  fpCommand->SetDefaultValue (6);
  fpCommand->SetRange("number >= 3");
}

G4VisCommandSetArrow3DLineSegmentsPerCircle::~G4VisCommandSetArrow3DLineSegmentsPerCircle ()
{
  delete fpCommand;
}

G4String G4VisCommandSetArrow3DLineSegmentsPerCircle::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetArrow3DLineSegmentsPerCircle::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  fCurrentArrow3DLineSegmentsPerCircle = fpCommand->GetNewIntValue(newValue);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout <<
    "Number of line segments per circle for drawing 3D arrows for future"
    "\n  \"/vis/scene/add/\" commands has been set to "
	   << fCurrentArrow3DLineSegmentsPerCircle
	   << G4endl;
  }
}

////////////// /vis/set/colour ////////////////////////////////////

G4VisCommandSetColour::G4VisCommandSetColour ()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/set/colour", this);
  fpCommand->SetGuidance
  ("Defines colour and opacity for future \"/vis/scene/add/\" commands.");
  fpCommand->SetGuidance
  ("(Except \"/vis/scene/add/text\" commands - see \"/vis/set/textColour\".)");
  fpCommand->SetGuidance(ConvertToColourGuidance());
  fpCommand->SetGuidance("Default: white and opaque.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("red", 's', omitable = true);
  parameter->SetGuidance
  ("Red component or a string, e.g., \"cyan\" (green and blue parameters are ignored).");
  parameter->SetDefaultValue ("1.");
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("alpha", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  parameter->SetGuidance ("Opacity");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSetColour::~G4VisCommandSetColour ()
{
  delete fpCommand;
}

G4String G4VisCommandSetColour::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetColour::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String redOrString;
  G4double green, blue, opacity;
  std::istringstream iss(newValue);
  iss >> redOrString >> green >> blue >> opacity;

  ConvertToColour(fCurrentColour, redOrString, green, blue, opacity);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout <<
    "Colour for future \"/vis/scene/add/\" commands has been set to "
    << fCurrentColour <<
    ".\n(Except \"/vis/scene/add/text\" commands - use \"/vis/set/textColour\".)"
    << G4endl;
  }
}

////////////// /vis/set/extentForField ////////////////////////////////////

G4VisCommandSetExtentForField::G4VisCommandSetExtentForField ()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/set/extentForField", this);
  fpCommand->SetGuidance
  ("Sets an extent for future \"/vis/scene/add/*Field\" commands.");
  fpCommand->SetGuidance
  ("The default is a null extent, which is interpreted by the commands as the"
   "\nextent of the whole scene.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("xmin", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("xmax", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("ymin", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("ymax", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("zmin", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("zmax", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSetExtentForField::~G4VisCommandSetExtentForField ()
{
  delete fpCommand;
}

G4String G4VisCommandSetExtentForField::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetExtentForField::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4double xmin, xmax, ymin, ymax, zmin, zmax;
  G4String unitString;
  std::istringstream iss(newValue);
  iss >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax >> unitString;
  G4double unit = G4UIcommand::ValueOf(unitString);
  xmin *= unit; xmax *= unit;
  ymin *= unit; ymax *= unit;
  zmin *= unit; zmax *= unit;

  fCurrentExtentForField = G4VisExtent(xmin,xmax,ymin,ymax,zmin,zmax);
  fCurrrentPVFindingsForField.clear();

  if (verbosity >= G4VisManager::confirmations) {
    G4cout <<
    "Extent for future \"/vis/scene/add/*Field\" commands has been set to "
    << fCurrentExtentForField
    << "\nVolume for field has been cleared."
    << G4endl;
  }
}

////////////// /vis/set/lineWidth ////////////////////////////////////

G4VisCommandSetLineWidth::G4VisCommandSetLineWidth ()
{
  G4bool omitable;
  fpCommand = new G4UIcmdWithADouble("/vis/set/lineWidth", this);
  fpCommand->SetGuidance
  ("Defines line width for future \"/vis/scene/add/\" commands.");
  fpCommand->SetParameterName ("lineWidth", omitable = true);
  fpCommand->SetDefaultValue (1.);
  fpCommand->SetRange("lineWidth >= 1.");
}

G4VisCommandSetLineWidth::~G4VisCommandSetLineWidth ()
{
  delete fpCommand;
}

G4String G4VisCommandSetLineWidth::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetLineWidth::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  fCurrentLineWidth = fpCommand->GetNewDoubleValue(newValue);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout <<
    "Line width for future \"/vis/scene/add/\" commands has been set to "
	   << fCurrentLineWidth
	   << G4endl;
  }
}

////////////// /vis/set/textColour ////////////////////////////////////

G4VisCommandSetTextColour::G4VisCommandSetTextColour ()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/set/textColour", this);
  fpCommand->SetGuidance
    ("Defines colour and opacity for future \"/vis/scene/add/text\" commands.");
  fpCommand->SetGuidance(ConvertToColourGuidance());
  fpCommand->SetGuidance("Default: blue and opaque.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("red", 's', omitable = true);
  parameter->SetGuidance
    ("Red component or a string, e.g., \"cyan\" (green and blue parameters are ignored).");
  parameter->SetDefaultValue ("0.");
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("alpha", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  parameter->SetGuidance ("Opacity");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSetTextColour::~G4VisCommandSetTextColour ()
{
  delete fpCommand;
}

G4String G4VisCommandSetTextColour::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetTextColour::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String redOrString;
  G4double green, blue, opacity;
  std::istringstream iss(newValue);
  iss >> redOrString >> green >> blue >> opacity;

  ConvertToColour(fCurrentTextColour, redOrString, green, blue, opacity);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout <<
      "Colour for future \"/vis/scene/add/text\" commands has been set to "
	   << fCurrentTextColour << '.'
	   << G4endl;
  }
}

////////////// /vis/set/textLayout ////////////////////////////////////

G4VisCommandSetTextLayout::G4VisCommandSetTextLayout ()
{
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/set/textLayout", this);
  fpCommand->SetGuidance
    ("Defines layout future \"/vis/scene/add/text\" commands.");
  fpCommand->SetGuidance
    ("\"left\" (default) for left justification to provided coordinate.");
  fpCommand->SetGuidance
    ("\"centre\" or \"center\" for text centered on provided coordinate.");
  fpCommand->SetGuidance
    ("\"right\" for right justification to provided coordinate.");
  fpCommand->SetGuidance("Default: left.");
  fpCommand->SetParameterName("layout", omitable = true);
  fpCommand->SetCandidates ("left centre center right");
  fpCommand->SetDefaultValue ("left");
}

G4VisCommandSetTextLayout::~G4VisCommandSetTextLayout ()
{
  delete fpCommand;
}

G4String G4VisCommandSetTextLayout::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetTextLayout::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4Text::Layout layout = G4Text::left;
  if (newValue == "left") layout = G4Text::left;
  else if (newValue == "centre" || newValue == "center")
    layout = G4Text::centre;
  else if (newValue == "right") layout = G4Text::right;

  fCurrentTextLayout = layout;

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Text layout (for future \"text\" commands) has been set to \""
	   << fCurrentTextLayout << "\"."
	   << G4endl;
  }
}

////////////// /vis/set/textSize ////////////////////////////////////

G4VisCommandSetTextSize::G4VisCommandSetTextSize ()
{
  G4bool omitable;
  fpCommand = new G4UIcmdWithADouble("/vis/set/textSize", this);
  fpCommand->SetGuidance
  ("Defines text size (pixels) for future \"/vis/scene/add/\" commands.");
  fpCommand->SetParameterName ("textSize", omitable = true);
  fpCommand->SetDefaultValue (12.);  // pixels
  fpCommand->SetRange("textSize >= 8.");
}

G4VisCommandSetTextSize::~G4VisCommandSetTextSize ()
{
  delete fpCommand;
}

G4String G4VisCommandSetTextSize::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetTextSize::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  fCurrentTextSize = fpCommand->GetNewDoubleValue(newValue);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout <<
    "Text size for future \"/vis/scene/add/\" commands has been set to "
    << fCurrentTextSize
    << G4endl;
  }
}

////////////// /vis/set/touchable ////////////////////////////////////

G4VisCommandSetTouchable::G4VisCommandSetTouchable ()
{
  G4bool omitable;
  G4UIparameter* parameter;
  fpCommand = new G4UIcommand("/vis/set/touchable", this);
  fpCommand->SetGuidance
  ("Defines touchable for future \"/vis/touchable/set/\" commands.");
  fpCommand->SetGuidance
  ("Please provide a list of space-separated physical volume names and"
   "\ncopy number pairs starting at the world volume, e.g:"
   "\n  /vis/set/touchable World 0 Envelope 0 Shape1 0"
   "\n(To get list of touchables, use \"/vis/drawTree\")"
   "\n(To save, use \"/vis/viewer/save\")");
  parameter = new G4UIparameter ("list", 's', omitable = true);
  parameter->SetGuidance
  ("List of physical volume names and copy number pairs");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSetTouchable::~G4VisCommandSetTouchable ()
{
  delete fpCommand;
}

G4String G4VisCommandSetTouchable::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetTouchable::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  if (newValue.isNull()) {
    fCurrentTouchableProperties = G4PhysicalVolumeModel::TouchableProperties();
    if (verbosity >= G4VisManager::confirmations) {
      G4cout <<
      "Current touchable reset to: " << fCurrentTouchableProperties.fTouchablePath
      << G4endl;
    }
    return;
  }
  
  G4ModelingParameters::PVNameCopyNoPath currentTouchablePath;

  // Algorithm from Josuttis p.476.
  G4String::size_type iBegin, iEnd;
  iBegin = newValue.find_first_not_of(' ');
  while (iBegin != G4String::npos) {
    iEnd = newValue.find_first_of(' ',iBegin);
    if (iEnd == G4String::npos) {
      iEnd = newValue.length();
    }
    G4String name(newValue.substr(iBegin,iEnd-iBegin));
    iBegin = newValue.find_first_not_of(' ',iEnd);
    if (iBegin == G4String::npos) {
      if (verbosity >= G4VisManager::warnings) {
        G4cout <<
        "WARNING: G4VisCommandSetTouchable::SetNewValue"
        "\n  A pair not found.  (Did you have an even number of parameters?)"
        "\n  Command ignored."
        << G4endl;
        return;
      }
    }
    iEnd = newValue.find_first_of(' ',iBegin);
    if (iEnd == G4String::npos) {
      iEnd = newValue.length();
    }
    G4int copyNo;
    std::istringstream iss(newValue.substr(iBegin,iEnd-iBegin));
    if (!(iss >> copyNo)) {
      if (verbosity >= G4VisManager::warnings) {
        G4cout <<
        "WARNING: G4VisCommandSetTouchable::SetNewValue"
        "\n  Error reading copy number - it was not numeric?"
        "\n  Command ignored."
        << G4endl;
        return;
      }
    }
    currentTouchablePath.push_back
    (G4ModelingParameters::PVNameCopyNo(name,copyNo));
    iBegin = newValue.find_first_not_of(' ',iEnd);
  }

  // Check validity
  G4bool successful = false;
  G4TransportationManager* transportationManager =
  G4TransportationManager::GetTransportationManager ();
  size_t nWorlds = transportationManager->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator iterWorld =
  transportationManager->GetWorldsIterator();
  for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
    G4PhysicalVolumeModel pvModel (*iterWorld);  // Unlimited depth.
    G4ModelingParameters mp;  // Default - no culling.
    pvModel.SetModelingParameters (&mp);
    G4TouchablePropertiesScene scene (&pvModel,currentTouchablePath);
    pvModel.DescribeYourselfTo (scene);  // Initiate geometry tree traversal.
    if (scene.GetFoundTouchableProperties().fpTouchablePV) {
      successful = true;
      fCurrentTouchableProperties = scene.GetFoundTouchableProperties();
      break;  // Found, so no need to scan more worlds.
    }
  }

  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout <<
      "Current touchable: " << fCurrentTouchableProperties.fTouchablePath
      << G4endl;
      return;
    }
  } else {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
      "WARNING: G4VisCommandSetTouchable::SetNewValue"
      "\n  Touchable not found."
      << G4endl;
      return;
    }
  }
}

////////////// /vis/set/volumeForField ////////////////////////////////////

G4VisCommandSetVolumeForField::G4VisCommandSetVolumeForField ()
{
  G4bool omitable;
  G4UIparameter* parameter;
  fpCommand = new G4UIcommand("/vis/set/volumeForField", this);
  fpCommand->SetGuidance
  ("Sets a volume for \"/vis/scene/add/*Field\" commands.");
  fpCommand->SetGuidance
  ("Takes a volume name or a /regular expression/ -- see guidance for"
   "\n\"/vis/drawVolume\"");
  parameter = new G4UIparameter ("physical-volume-name", 's', omitable = true);
  parameter -> SetDefaultValue ("none");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("copy-no", 'i', omitable = true);
  parameter -> SetGuidance ("If negative, matches any copy no.");
  parameter -> SetDefaultValue (-1);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("draw", 'b', omitable = true);
  parameter -> SetGuidance ("If true, draw extent of found volumes.");
  parameter -> SetDefaultValue (false);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSetVolumeForField::~G4VisCommandSetVolumeForField ()
{
  delete fpCommand;
}

G4String G4VisCommandSetVolumeForField::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetVolumeForField::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String name, drawString;
  G4int copyNo;
  std::istringstream is (newValue);
  is >> name >> copyNo >> drawString;
  G4bool draw = G4UIcmdWithABool::ConvertToBool(drawString);

  if (name == "none") {
    fCurrrentPVFindingsForField.clear();
    fCurrentExtentForField = G4VisExtent();
    if (verbosity >= G4VisManager::warnings) {
      G4cout << "Volume for field cleared" << G4endl;
    }
    return;
  }

  G4TransportationManager* transportationManager =
  G4TransportationManager::GetTransportationManager ();
  size_t nWorlds = transportationManager->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator iterWorld =
  transportationManager->GetWorldsIterator();
  fCurrrentPVFindingsForField.clear();
  G4BoundingExtentScene extentScene;
  for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
    G4PhysicalVolumeModel searchModel (*iterWorld);  // Unlimited depth.
    G4ModelingParameters mp;  // Default - no culling.
    searchModel.SetModelingParameters (&mp);
    G4PhysicalVolumesSearchScene searchScene (&searchModel, name, copyNo);
    searchModel.DescribeYourselfTo (searchScene);  // Initiate search.
    for (const auto& findings: searchScene.GetFindings()) {
      fCurrrentPVFindingsForField.push_back(findings);
      G4VisExtent extent = findings.fpFoundPV->GetLogicalVolume()->GetSolid()->GetExtent();
      extent.Transform(findings.fFoundObjectTransformation);
      extentScene.AccrueBoundingExtent(extent);
    }
  }

  if (fCurrrentPVFindingsForField.empty()) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Volume \"" << name << "\"";
      if (copyNo >= 0) {
        G4cerr << ", copy no. " << copyNo << ",";
      }
      G4cerr << " not found." << G4endl;
    }
    return;
  }

  fCurrentExtentForField = extentScene.GetExtent();

  if (draw) DrawExtent(fCurrentExtentForField);

  if (verbosity >= G4VisManager::confirmations) {
    for (const auto& findings: fCurrrentPVFindingsForField) {
      G4cout
      << "\"" << findings.fpFoundPV->GetName()
      << "\", copy no. " << findings.fFoundPVCopyNo
      << ", found\nin searched volume \""
      << findings.fpSearchPV->GetName()
      << "\" at depth " << findings.fFoundDepth
      << ",\nbase path: \"" << findings.fFoundBasePVPath
      << "\",\nand has been set as volume for field."
      << G4endl;
    }
  }
}
