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
// /vis/scene/add commands - John Allison  9th August 1998

#include "G4VisCommandsSceneAdd.hh"

#include "G4TransportationManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4ModelingParameters.hh"
#include "G4HitsModel.hh"
#include "G4DigiModel.hh"
#include "G4GPSModel.hh"
#include "G4ElectricFieldModel.hh"
#include "G4MagneticFieldModel.hh"
#include "G4PSHitsModel.hh"
#include "G4TrajectoriesModel.hh"
#include "G4TextModel.hh"
#include "G4ArrowModel.hh"
#include "G4AxesModel.hh"
#include "G4PlotterModel.hh"
#include "G4PhysicalVolumesSearchScene.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ApplicationState.hh"
#include "G4VUserVisAction.hh"
#include "G4CallbackModel.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Polyhedron.hh"
#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4Tokenizer.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4StateManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "G4RichTrajectory.hh"
#include "G4RichTrajectoryPoint.hh"
#include "G4SmoothTrajectory.hh"
#include "G4SmoothTrajectoryPoint.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"
#include "G4Polyline.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSourceData.hh"
#include "G4PlotterManager.hh"

#include <sstream>

#define G4warn G4cout

////////////// /vis/scene/add/arrow ///////////////////////////////////////

G4VisCommandSceneAddArrow::G4VisCommandSceneAddArrow () {
  fpCommand = new G4UIcommand("/vis/scene/add/arrow", this);
  fpCommand -> SetGuidance ("Adds arrow to current scene.");
  G4bool omitable;
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("x1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("z1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("x2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("z2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter    (parameter);
}

G4VisCommandSceneAddArrow::~G4VisCommandSceneAddArrow () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddArrow::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddArrow::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String unitString;
  G4double x1, y1, z1, x2, y2, z2;
  std::istringstream is(newValue);
  is >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> unitString;
  G4double unit = G4UIcommand::ValueOf(unitString);
  x1 *= unit; y1 *= unit; z1 *= unit;
  x2 *= unit; y2 *= unit; z2 *= unit;

  // Consult scene for arrow width.
  const G4VisExtent& sceneExtent = pScene->GetExtent();
  G4double arrowWidth =
    0.005 * fCurrentLineWidth * sceneExtent.GetExtentRadius();

  G4VModel* model = new G4ArrowModel
    (x1, y1, z1, x2, y2, z2,
     arrowWidth, fCurrentColour, newValue,
     fCurrentArrow3DLineSegmentsPerCircle);

  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Arrow has been added to scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/add/arrow2D ///////////////////////////////////////

G4VisCommandSceneAddArrow2D::G4VisCommandSceneAddArrow2D () {
  fpCommand = new G4UIcommand("/vis/scene/add/arrow2D", this);
  fpCommand -> SetGuidance ("Adds 2D arrow to current scene.");
  G4bool omitable;
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("x1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("x2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddArrow2D::~G4VisCommandSceneAddArrow2D () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddArrow2D::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddArrow2D::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4double x1, y1, x2, y2;
  std::istringstream is(newValue);
  is >> x1 >> y1 >> x2 >> y2;

  Arrow2D* arrow2D = new Arrow2D
    (x1, y1, x2, y2, fCurrentLineWidth, fCurrentColour);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddArrow2D::Arrow2D>(arrow2D);
  model->SetType("Arrow2D");
  model->SetGlobalTag("Arrow2D");
  model->SetGlobalDescription("Arrow2D: " + newValue);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "A 2D arrow has been added to scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

G4VisCommandSceneAddArrow2D::Arrow2D::Arrow2D
(G4double x1, G4double y1,
 G4double x2, G4double y2,
 G4double width, const G4Colour& colour):
  fWidth(width), fColour(colour)
{
  fShaftPolyline.push_back(G4Point3D(x1,y1,0));
  fShaftPolyline.push_back(G4Point3D(x2,y2,0));
  G4Vector3D arrowDirection = G4Vector3D(x2-x1,y2-y1,0).unit();
  G4Vector3D arrowPointLeftDirection(arrowDirection);
  arrowPointLeftDirection.rotateZ(150.*deg);
  G4Vector3D arrowPointRightDirection(arrowDirection);
  arrowPointRightDirection.rotateZ(-150.*deg);
  fHeadPolyline.push_back(G4Point3D(x2,y2,0)+0.04*arrowPointLeftDirection);
  fHeadPolyline.push_back(G4Point3D(x2,y2,0));
  fHeadPolyline.push_back(G4Point3D(x2,y2,0)+0.04*arrowPointRightDirection);
  G4VisAttributes va;
  va.SetLineWidth(fWidth);
  va.SetColour(fColour);
  fShaftPolyline.SetVisAttributes(va);
  fHeadPolyline.SetVisAttributes(va);
}

void G4VisCommandSceneAddArrow2D::Arrow2D::operator()
  (G4VGraphicsScene& sceneHandler, const G4ModelingParameters*)
{
  sceneHandler.BeginPrimitives2D();
  sceneHandler.AddPrimitive(fShaftPolyline);
  sceneHandler.AddPrimitive(fHeadPolyline);
  sceneHandler.EndPrimitives2D();
}

////////////// /vis/scene/add/axes //////////////////////////////////

G4VisCommandSceneAddAxes::G4VisCommandSceneAddAxes () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/axes", this);
  fpCommand -> SetGuidance ("Add axes.");
  fpCommand -> SetGuidance
  ("Draws axes at (x0, y0, z0) of given length and colour.");
  fpCommand -> SetGuidance
  ("If \"colour-string\" is \"auto\", x, y and z will be red, green and blue"
   "\n  respectively.  Otherwise it can be one of the pre-defined text-specified"
   "\n  colours - see information printed by the vis manager at start-up or"
   "\n  use \"/vis/list\".");
  fpCommand -> SetGuidance
  ("If \"length\" is negative, it is set to about 25% of scene extent.");
  fpCommand -> SetGuidance
  ("If \"showtext\" is false, annotations are suppressed.");
  G4UIparameter* parameter;
  parameter =  new G4UIparameter ("x0", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("y0", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("z0", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("length", 'd', omitable = true);
  parameter->SetDefaultValue (-1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("colour-string", 's', omitable = true);
  parameter->SetDefaultValue  ("auto");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("showtext", 'b', omitable = true);
  parameter->SetDefaultValue  ("true");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSceneAddAxes::~G4VisCommandSceneAddAxes () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddAxes::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddAxes::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  } else {
    if (pScene->GetExtent().GetExtentRadius() <= 0.) {
      if (verbosity >= G4VisManager::errors) {
        G4warn
  << "ERROR: Scene has no extent. Add volumes or use \"/vis/scene/add/extent\"."
        << G4endl;
      }
      return;
    }
  }

  G4String unitString, colourString, showTextString;
  G4double x0, y0, z0, length;
  std::istringstream is (newValue);
  is >> x0 >> y0 >> z0 >> length >> unitString
  >> colourString >> showTextString;
  G4bool showText = G4UIcommand::ConvertToBool(showTextString);


  G4double unit = G4UIcommand::ValueOf(unitString);
  x0 *= unit; y0 *= unit; z0 *= unit;
  const G4VisExtent& sceneExtent = pScene->GetExtent();  // Existing extent.
  if (length < 0.) {
    const G4double lengthMax = 0.5 * sceneExtent.GetExtentRadius();
    const G4double intLog10Length = std::floor(std::log10(lengthMax));
    length = std::pow(10,intLog10Length);
    if (5.*length < lengthMax) length *= 5.;
    else if (2.*length < lengthMax) length *= 2.;
  } else {
    length *= unit;
  }

  // Consult scene for arrow width...
  G4double arrowWidth =
    0.05 * fCurrentLineWidth * sceneExtent.GetExtentRadius();
  // ...but limit it to length/30.
  if (arrowWidth > length/30.) arrowWidth = length/30.;

  G4VModel* model = new G4AxesModel
    (x0, y0, z0, length, arrowWidth, colourString, newValue,
     showText, fCurrentTextSize);

  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  const G4String& currentSceneName = pScene -> GetName ();
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Axes of length " << G4BestUnit(length,"Length")
      << "have been added to scene \"" << currentSceneName << "\"."
      << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/add/date ///////////////////////////////////////

G4VisCommandSceneAddDate::G4VisCommandSceneAddDate () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/date", this);
  fpCommand -> SetGuidance ("Adds date to current scene.");
  fpCommand -> SetGuidance
  ("If \"date\"is omitted, the current date and time is drawn."
   "\nOtherwise, the string, including the rest of the line, is drawn.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("size", 'i', omitable = true);
  parameter -> SetGuidance ("Screen size of text in pixels.");
  parameter -> SetDefaultValue (18);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("x-position", 'd', omitable = true);
  parameter -> SetGuidance ("x screen position in range -1 < x < 1.");
  parameter -> SetDefaultValue (0.95);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y-position", 'd', omitable = true);
  parameter -> SetGuidance ("y screen position in range -1 < y < 1.");
  parameter -> SetDefaultValue (0.9);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("layout", 's', omitable = true);
  parameter -> SetGuidance ("Layout, i.e., adjustment: left|centre|right.");
  parameter -> SetDefaultValue ("right");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("date", 's', omitable = true);
  parameter -> SetDefaultValue ("-");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddDate::~G4VisCommandSceneAddDate () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddDate::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddDate::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int size;
  G4double x, y;
  G4String layoutString, dateString;
  std::istringstream is(newValue);
  is >> size >> x >> y >> layoutString >> dateString;
  // Read rest of line, if any.
  const size_t NREMAINDER = 100;
  char remainder[NREMAINDER];
  remainder[0]='\0';  // In case there is nothing remaining.
  is.getline(remainder, NREMAINDER);
  dateString += remainder;
  G4Text::Layout layout = G4Text::right;
  if (layoutString[0] == 'l') layout = G4Text::left;
  else if (layoutString[0] == 'c') layout = G4Text::centre;
  else if (layoutString[0] == 'r') layout = G4Text::right;

  Date* date = new Date(fpVisManager, size, x, y, layout, dateString);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddDate::Date>(date);
  model->SetType("Date");
  model->SetGlobalTag("Date");
  model->SetGlobalDescription("Date: " + newValue);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Date has been added to scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

void G4VisCommandSceneAddDate::Date::operator()
  (G4VGraphicsScene& sceneHandler, const G4ModelingParameters*)
{
  G4String time;
  if (fDate == "-") {
    time = fTimer.GetClockTime();
  } else {
    time = fDate;
  }
  // Check for \n, starting from back, and erase.
  std::string::size_type i = time.rfind('\n');
  if (i != std::string::npos) time.erase(i);
  G4Text text(time, G4Point3D(fX, fY, 0.));
  text.SetScreenSize(fSize);
  text.SetLayout(fLayout);
  G4VisAttributes textAtts(G4Colour(0.,1.,1));
  text.SetVisAttributes(textAtts);
  sceneHandler.BeginPrimitives2D();
  sceneHandler.AddPrimitive(text);
  sceneHandler.EndPrimitives2D();
}

////////////// /vis/scene/add/digis ///////////////////////////////////////

G4VisCommandSceneAddDigis::G4VisCommandSceneAddDigis () {
  fpCommand = new G4UIcmdWithoutParameter ("/vis/scene/add/digis", this);
  fpCommand -> SetGuidance ("Adds digis to current scene.");
  fpCommand -> SetGuidance
    ("Digis are drawn at end of event when the scene in which"
     "\nthey are added is current.");
}

G4VisCommandSceneAddDigis::~G4VisCommandSceneAddDigis () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddDigis::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddDigis::SetNewValue (G4UIcommand*, G4String) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4VModel* model = new G4DigiModel;
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddEndOfEventModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Digis, if any, will be drawn at end of run in scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/add/electricField ///////////////////////////////////////

G4VisCommandSceneAddElectricField::G4VisCommandSceneAddElectricField () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/electricField", this);
  fpCommand -> SetGuidance
  ("Adds electric field representation to current scene.");
  fpCommand -> SetGuidance
  ("The first parameter is no. of data points per half extent.  So, possibly, at"
   "\nmaximum, the number of data points sampled is (2*n+1)^3, which can grow"
   "\nlarge--be warned!"
   "\nThe default value is 10, i.e., a 21x21x21 array, i.e., 9,261 sampling points."
   "\nThat may swamp your view, but usually, a field is limited to a small part of"
   "\nthe extent, so it's not a problem. But if it is, here are some of the things"
   "\nyou can do:"
   "\n- reduce the number of data points per half extent (first parameter);"
   "\n- specify \"lightArrow\" (second parameter);"
   "\n- restrict the region sampled with \"/vis/set/extentForField\";"
   "\n- restrict the drawing to a specific volume with"
   "\n    \"/vis/set/volumeForField\" or \"/vis/touchable/volumeForField\"."
   "\nNote: you might have to deactivate existing field models with"
   "\n  \"/vis/scene/activateModel Field false\" and re-issue"
   "\n  \"/vis/scene/add/...Field\" command again.");
  fpCommand -> SetGuidance
  ("In the arrow representation, the length of the arrow is proportional"
   "\nto the magnitude of the field and the colour is mapped onto the range"
   "\nas a fraction of the maximum magnitude: 0->0.5->1 is red->green->blue.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("nDataPointsPerHalfExtent", 'i', omitable = true);
  parameter -> SetDefaultValue (10);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("representation", 's', omitable = true);
  parameter -> SetParameterCandidates("fullArrow lightArrow");
  parameter -> SetDefaultValue ("fullArrow");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddElectricField::~G4VisCommandSceneAddElectricField () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddElectricField::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddElectricField::SetNewValue
(G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<  "ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int nDataPointsPerHalfExtent;
  G4String representation;
  std::istringstream iss(newValue);
  iss >> nDataPointsPerHalfExtent >> representation;
  G4ElectricFieldModel::Representation
  modelRepresentation = G4ElectricFieldModel::fullArrow;
  if (representation == "lightArrow") {
    modelRepresentation = G4ElectricFieldModel::lightArrow;
  }
  G4VModel* model;
  model = new G4ElectricFieldModel
  (nDataPointsPerHalfExtent,modelRepresentation,
   fCurrentArrow3DLineSegmentsPerCircle,
   fCurrentExtentForField,
   fCurrrentPVFindingsForField);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout
      << "Electric field, if any, will be drawn in scene \""
      << currentSceneName
      << "\"\n  with "
      << nDataPointsPerHalfExtent
      << " data points per half extent and with representation \""
      << representation
      << '\"'
      << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/add/eventID ///////////////////////////////////////

G4VisCommandSceneAddEventID::G4VisCommandSceneAddEventID () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/eventID", this);
  fpCommand -> SetGuidance ("Adds eventID to current scene.");
  fpCommand -> SetGuidance
    ("Run and event numbers are drawn at end of event or run when"
     "\n the scene in which they are added is current.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("size", 'i', omitable = true);
  parameter -> SetGuidance ("Screen size of text in pixels.");
  parameter -> SetDefaultValue (18);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("x-position", 'd', omitable = true);
  parameter -> SetGuidance ("x screen position in range -1 < x < 1.");
  parameter -> SetDefaultValue (-0.95);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y-position", 'd', omitable = true);
  parameter -> SetGuidance ("y screen position in range -1 < y < 1.");
  parameter -> SetDefaultValue (0.9);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("layout", 's', omitable = true);
  parameter -> SetGuidance ("Layout, i.e., adjustment: left|centre|right.");
  parameter -> SetDefaultValue ("left");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddEventID::~G4VisCommandSceneAddEventID () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddEventID::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddEventID::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int size;
  G4double x, y;
  G4String layoutString;
  std::istringstream is(newValue);
  is >> size >> x >> y >> layoutString;

  G4Text::Layout layout = G4Text::right;
  if (layoutString[0] == 'l') layout = G4Text::left;
  else if (layoutString[0] == 'c') layout = G4Text::centre;
  else if (layoutString[0] == 'r') layout = G4Text::right;

  // For End of Event (only for reviewing kept events one by one)
  EventID* eoeEventID
  = new EventID(forEndOfEvent, fpVisManager, size, x, y, layout);
  G4VModel* eoeModel =
    new G4CallbackModel<G4VisCommandSceneAddEventID::EventID>(eoeEventID);
  eoeModel->SetType("EoEEventID");
  eoeModel->SetGlobalTag("EoEEventID");
  eoeModel->SetGlobalDescription("EoEEventID: " + newValue);
  G4bool successfulEoE = pScene -> AddEndOfEventModel (eoeModel, warn);

  // For End of Run
  EventID* eorEventID
  = new EventID(forEndOfRun, fpVisManager, size, x, y, layout);
  G4VModel* eorModel =
  new G4CallbackModel<G4VisCommandSceneAddEventID::EventID>(eorEventID);
  eorModel->SetType("EoREventID");
  eorModel->SetGlobalTag("EoREventID");
  eorModel->SetGlobalDescription("EoREventID: " + newValue);
  G4bool successfulEoR = pScene -> AddEndOfRunModel (eorModel, warn);

  if (successfulEoE && successfulEoR) {
    if (verbosity >= G4VisManager::confirmations) {
      const G4String& currentSceneName = pScene -> GetName ();
      G4cout << "EventID has been added to scene \""
      << currentSceneName << "\"."
      << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

void G4VisCommandSceneAddEventID::EventID::operator()
(G4VGraphicsScene& sceneHandler, const G4ModelingParameters* mp)
{
  G4RunManager* runManager = G4RunManagerFactory::GetMasterRunManager();
  if(!runManager)
    return;

  const G4Run* currentRun = runManager->GetCurrentRun();
  if (!currentRun) return;

  const G4int currentRunID = currentRun->GetRunID();

  std::ostringstream oss;
  switch (fForWhat) {
    case forEndOfEvent:
    {
      // Only use if reviewing kept events
      if (!fpVisManager->GetReviewingKeptEvents()) return;
      const G4Event* currentEvent = mp->GetEvent();
      if (!currentEvent) return;
      G4int eventID = currentEvent->GetEventID();
      oss << "Run " << currentRunID << " Event " << eventID;
      break;
    }
    case forEndOfRun:
    {
      // Only use if NOT reviewing kept events
      if (fpVisManager->GetReviewingKeptEvents()) return;
      const G4int nEvents = currentRun->GetNumberOfEventToBeProcessed();
      const auto* events = currentRun->GetEventVector();
      size_t nKeptEvents = events? events->size(): 0;
      oss << "Run " << currentRunID << " (" << nEvents << " event";
      if (nEvents != 1) oss << 's';
      oss << ", " << nKeptEvents << " kept)";
      break;
    }
    default:
      return;
  }

  G4Text text(oss.str(), G4Point3D(fX, fY, 0.));
  text.SetScreenSize(fSize);
  text.SetLayout(fLayout);
  G4VisAttributes textAtts(G4Colour(0.,1.,1));
  text.SetVisAttributes(textAtts);
  sceneHandler.BeginPrimitives2D();
  sceneHandler.AddPrimitive(text);
  sceneHandler.EndPrimitives2D();
}

////////////// /vis/scene/add/extent ///////////////////////////////////////

G4VisCommandSceneAddExtent::G4VisCommandSceneAddExtent () {
  fpCommand = new G4UIcommand("/vis/scene/add/extent", this);
  fpCommand -> SetGuidance
  ("Adds a dummy model with given extent to the current scene."
   "\nRequires the limits: xmin, xmax, ymin, ymax, zmin, zmax unit."
   "\nThis can be used to provide an extent to the scene even if"
   "\nno other models with extent are available. For example,"
   "\neven if there is no geometry.  In that case, for example:"
   "\n  /vis/open OGL"
   "\n  /vis/scene/create"
   "\n  /vis/scene/add/extent -300 300 -300 300 -300 300 cm"
   "\n  /vis/sceneHandler/attach");
  G4bool omitable;
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("xmin", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("xmax", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("ymin", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("ymax", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("zmin", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("zmax", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter -> SetDefaultValue ("m");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddExtent::~G4VisCommandSceneAddExtent () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddExtent::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddExtent::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4double xmin, xmax, ymin, ymax, zmin, zmax;
  G4String unitString;
  std::istringstream is(newValue);
  is >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax >> unitString;
  G4double unit = G4UIcommand::ValueOf(unitString);
  xmin *= unit; xmax *= unit;
  ymin *= unit; ymax *= unit;
  zmin *= unit; zmax *= unit;

  G4VisExtent visExtent(xmin, xmax, ymin, ymax, zmin, zmax);
  Extent* extent = new Extent(xmin, xmax, ymin, ymax, zmin, zmax);
  G4VModel* model =
  new G4CallbackModel<G4VisCommandSceneAddExtent::Extent>(extent);
  model->SetType("Extent");
  model->SetGlobalTag("Extent");
  model->SetGlobalDescription("Extent: " + newValue);
  model->SetExtent(visExtent);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "A benign model with extent "
      << visExtent
      << " has been added to scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

G4VisCommandSceneAddExtent::Extent::Extent
(G4double xmin, G4double xmax,
 G4double ymin, G4double ymax,
 G4double zmin, G4double zmax):
fExtent(xmin,xmax,ymin,ymax,zmin,zmax)
{}

void G4VisCommandSceneAddExtent::Extent::operator()
(G4VGraphicsScene&, const G4ModelingParameters*)
{}

////////////// /vis/scene/add/frame ///////////////////////////////////////

G4VisCommandSceneAddFrame::G4VisCommandSceneAddFrame () {
  fpCommand = new G4UIcommand("/vis/scene/add/frame", this);
  fpCommand -> SetGuidance ("Add frame to current scene.");
  G4bool omitable;
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("size", 'd', omitable = true);
  parameter -> SetGuidance ("Size of frame.  1 = full window.");
  parameter -> SetParameterRange ("size > 0 && size <=1");
  parameter -> SetDefaultValue (0.97);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddFrame::~G4VisCommandSceneAddFrame () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddFrame::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddFrame::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4double size;
  std::istringstream is(newValue);
  is >> size;

  Frame* frame = new Frame(size, fCurrentLineWidth, fCurrentColour);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddFrame::Frame>(frame);
  model->SetType("Frame");
  model->SetGlobalTag("Frame");
  model->SetGlobalDescription("Frame: " + newValue);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Frame has been added to scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

void G4VisCommandSceneAddFrame::Frame::operator()
  (G4VGraphicsScene& sceneHandler, const G4ModelingParameters*)
{
  G4Polyline frame;
  frame.push_back(G4Point3D( fSize,  fSize, 0.));
  frame.push_back(G4Point3D(-fSize,  fSize, 0.));
  frame.push_back(G4Point3D(-fSize, -fSize, 0.));
  frame.push_back(G4Point3D( fSize, -fSize, 0.));
  frame.push_back(G4Point3D( fSize,  fSize, 0.));
  G4VisAttributes va;
  va.SetLineWidth(fWidth);
  va.SetColour(fColour);
  frame.SetVisAttributes(va);
  sceneHandler.BeginPrimitives2D();
  sceneHandler.AddPrimitive(frame);
  sceneHandler.EndPrimitives2D();
}

////////////// /vis/scene/add/gps ///////////////////////////////////////

G4VisCommandSceneAddGPS::G4VisCommandSceneAddGPS () {
  G4bool omitable;
  G4UIparameter* parameter;
  fpCommand = new G4UIcommand ("/vis/scene/add/gps", this);
  fpCommand -> SetGuidance
  ("A representation of the source(s) of the General Particle Source"
   "\nwill be added to current scene and drawn, if applicable.");
  fpCommand->SetGuidance(ConvertToColourGuidance());
  fpCommand->SetGuidance("Default: red and transparent.");
  parameter = new G4UIparameter("red_or_string", 's', omitable = true);
  parameter -> SetDefaultValue ("1.");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("green", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("opacity", 'd', omitable = true);
  parameter -> SetDefaultValue (0.3);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddGPS::~G4VisCommandSceneAddGPS () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddGPS::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddGPS::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String redOrString;
  G4double green, blue, opacity;
  std::istringstream iss(newValue);
  iss >> redOrString >> green >> blue >> opacity;
  G4Colour colour(1.,0.,0.,0.3);  // Default red and transparent.
  ConvertToColour(colour, redOrString, green, blue, opacity);

  G4VModel* model = new G4GPSModel(colour);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout <<
      "A representation of the source(s) of the General Particle Source will be drawn"
      "\n  in colour " << colour << " for scene \""
      << currentSceneName << "\" if applicable."
      << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}
  
////////////// /vis/scene/add/hits ///////////////////////////////////////

G4VisCommandSceneAddHits::G4VisCommandSceneAddHits () {
  fpCommand = new G4UIcmdWithoutParameter ("/vis/scene/add/hits", this);
  fpCommand -> SetGuidance ("Adds hits to current scene.");
  fpCommand -> SetGuidance
    ("Hits are drawn at end of event when the scene in which"
     "\nthey are added is current.");
}

G4VisCommandSceneAddHits::~G4VisCommandSceneAddHits () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddHits::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddHits::SetNewValue (G4UIcommand*, G4String) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4VModel* model = new G4HitsModel;
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddEndOfEventModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Hits, if any, will be drawn at end of run in scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/add/line ///////////////////////////////////////

G4VisCommandSceneAddLine::G4VisCommandSceneAddLine () {
  fpCommand = new G4UIcommand("/vis/scene/add/line", this);
  fpCommand -> SetGuidance ("Adds line to current scene.");
  G4bool omitable;
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("x1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("z1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("x2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("z2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter    (parameter);
}

G4VisCommandSceneAddLine::~G4VisCommandSceneAddLine () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddLine::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddLine::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String unitString;
  G4double x1, y1, z1, x2, y2, z2;
  std::istringstream is(newValue);
  is >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> unitString;
  G4double unit = G4UIcommand::ValueOf(unitString);
  x1 *= unit; y1 *= unit; z1 *= unit;
  x2 *= unit; y2 *= unit; z2 *= unit;

  Line* line = new Line(x1, y1, z1, x2, y2, z2,
			fCurrentLineWidth, fCurrentColour);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddLine::Line>(line);
  model->SetType("Line");
  model->SetGlobalTag("Line");
  model->SetGlobalDescription("Line: " + newValue);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Line has been added to scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

G4VisCommandSceneAddLine::Line::Line
(G4double x1, G4double y1, G4double z1,
 G4double x2, G4double y2, G4double z2,
 G4double width, const G4Colour& colour):
  fWidth(width), fColour(colour)
{
  fPolyline.push_back(G4Point3D(x1,y1,z1));
  fPolyline.push_back(G4Point3D(x2,y2,z2));
  G4VisAttributes va;
  va.SetLineWidth(fWidth);
  va.SetColour(fColour);
  fPolyline.SetVisAttributes(va);
}

void G4VisCommandSceneAddLine::Line::operator()
  (G4VGraphicsScene& sceneHandler, const G4ModelingParameters*)
{
  sceneHandler.BeginPrimitives();
  sceneHandler.AddPrimitive(fPolyline);
  sceneHandler.EndPrimitives();
}

////////////// /vis/scene/add/line2D ///////////////////////////////////////

G4VisCommandSceneAddLine2D::G4VisCommandSceneAddLine2D () {
  fpCommand = new G4UIcommand("/vis/scene/add/line2D", this);
  fpCommand -> SetGuidance ("Adds 2D line to current scene.");
  G4bool omitable;
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("x1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y1", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("x2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y2", 'd', omitable = false);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddLine2D::~G4VisCommandSceneAddLine2D () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddLine2D::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddLine2D::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4double x1, y1, x2, y2;
  std::istringstream is(newValue);
  is >> x1 >> y1 >> x2 >> y2;

  Line2D* line2D = new Line2D
    (x1, y1, x2, y2, fCurrentLineWidth, fCurrentColour);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddLine2D::Line2D>(line2D);
  model->SetType("Line2D");
  model->SetGlobalTag("Line2D");
  model->SetGlobalDescription("Line2D: " + newValue);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "A 2D line has been added to scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

G4VisCommandSceneAddLine2D::Line2D::Line2D
(G4double x1, G4double y1,
 G4double x2, G4double y2,
 G4double width, const G4Colour& colour):
  fWidth(width), fColour(colour)
{
  fPolyline.push_back(G4Point3D(x1,y1,0));
  fPolyline.push_back(G4Point3D(x2,y2,0));
  G4VisAttributes va;
  va.SetLineWidth(fWidth);
  va.SetColour(fColour);
  fPolyline.SetVisAttributes(va);
}

void G4VisCommandSceneAddLine2D::Line2D::operator()
  (G4VGraphicsScene& sceneHandler, const G4ModelingParameters*)
{
  sceneHandler.BeginPrimitives2D();
  sceneHandler.AddPrimitive(fPolyline);
  sceneHandler.EndPrimitives2D();
}

////////////// /vis/scene/add/localAxes ///////////////////////////////////////

G4VisCommandSceneAddLocalAxes::G4VisCommandSceneAddLocalAxes () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/localAxes", this);
  fpCommand -> SetGuidance
  ("Adds local axes to physical volume(s).");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("physical-volume-name", 's', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("copy-no", 'i', omitable = true);
  parameter -> SetGuidance ("If negative, matches any copy no.");
  parameter -> SetDefaultValue (-1);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddLocalAxes::~G4VisCommandSceneAddLocalAxes () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddLocalAxes::GetCurrentValue (G4UIcommand*) {
  return "world 0 -1";
}

void G4VisCommandSceneAddLocalAxes::SetNewValue (G4UIcommand*,
						 G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn = verbosity >= G4VisManager::warnings;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String name;
  G4int copyNo;
  std::istringstream is (newValue);
  is >> name >> copyNo;

  std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVector;

  // Search all worlds...
  G4TransportationManager* transportationManager =
  G4TransportationManager::GetTransportationManager ();
  std::vector<G4VPhysicalVolume*>::iterator iterWorld =
  transportationManager->GetWorldsIterator();
  size_t nWorlds = transportationManager->GetNoWorlds();
  for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
    G4ModelingParameters mp;  // Default - no culling.
    G4PhysicalVolumeModel searchModel
    (*iterWorld,
     G4PhysicalVolumeModel::UNLIMITED,
     G4Transform3D(),
     &mp,
     true);  // Use full extent (avoids initial descent of geometry tree)
    G4PhysicalVolumesSearchScene searchScene
    (&searchModel, name, copyNo);
    searchModel.DescribeYourselfTo (searchScene);  // Initiate search.
    for (const auto& findings: searchScene.GetFindings()) {
      findingsVector.push_back(findings);
    }
  }

  G4int id = 0;  // To distinguish axes models by their global description
  for (const auto& findings: findingsVector) {

    // Create axes model based on size and transformation of found volume(s).
    const auto& extent = findings.fpFoundPV->GetLogicalVolume()->GetSolid()->GetExtent();
    const auto& transform = findings.fFoundObjectTransformation;

    const G4double lengthMax = extent.GetExtentRadius()/2.;
    const G4double intLog10LengthMax = std::floor(std::log10(lengthMax));
    G4double length = std::pow(10,intLog10LengthMax);
    if (5.*length < lengthMax) length *= 5.;
    else if (2.*length < lengthMax) length *= 2.;

    const auto& axesModel = new G4AxesModel(0.,0.,0.,length,transform);
    axesModel->SetGlobalTag("LocalAxesModel");
    std::ostringstream oss; oss
    << "Local Axes for " << findings.fpFoundPV->GetName()
    << ':' << findings.fFoundPVCopyNo << ':' << id++;
    axesModel->SetGlobalDescription(oss.str());
    // ...so add it to the scene.
    G4bool successful = pScene->AddRunDurationModel(axesModel,warn);
    if (successful) {
      if (verbosity >= G4VisManager::confirmations) {
	G4cout << "\"" << findings.fpFoundPV->GetName()
	<< "\", copy no. " << findings.fFoundPVCopyNo
	<< ",\n  found in searched volume \""
	<< findings.fpSearchPV->GetName()
	<< "\" at depth " << findings.fFoundDepth
	<< ",\n  base path: \"" << findings.fFoundBasePVPath
	<< "\".\n  Local axes have been added to scene \""
	<< pScene->GetName() << "\".";
	if (verbosity >= G4VisManager::parameters) {
	  G4cout << "  With extent " << extent
	  << "\n  at " << transform.getRotation()
	  << "  " << transform.getTranslation();
	}
	G4cout << G4endl;
      }
    } else {
      G4VisCommandsSceneAddUnsuccessful(verbosity);
    }
  }

  if (findingsVector.empty()) {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: Volume \"" << name << "\"";
      if (copyNo >= 0) {
	G4warn << ", copy no. " << copyNo << ",";
      }
      G4warn << " not found." << G4endl;
    }
    G4VisCommandsSceneAddUnsuccessful(verbosity);
    return;
  }

  CheckSceneAndNotifyHandlers(pScene);
}

////////////// /vis/scene/add/logicalVolume //////////////////////////////////

G4VisCommandSceneAddLogicalVolume::G4VisCommandSceneAddLogicalVolume () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/logicalVolume", this);
  fpCommand -> SetGuidance ("Adds a logical volume to the current scene,");
  fpCommand -> SetGuidance
  ("Shows boolean components (if any), voxels (if any), readout geometry"
   "\n  (if any), local axes and overlaps (if any), under control of the"
   "\n  appropriate flag."
   "\n  Note: voxels are not constructed until start of run -"
   "\n \"/run/beamOn\".  (For voxels without a run, \"/run/beamOn 0\".)");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("logical-volume-name", 's', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("depth-of-descent", 'i', omitable = true);
  parameter -> SetGuidance ("Depth of descent of geometry hierarchy.");
  parameter -> SetDefaultValue (1);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("booleans-flag", 'b', omitable = true);
  parameter -> SetDefaultValue (true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("voxels-flag", 'b', omitable = true);
  parameter -> SetDefaultValue (true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("readout-flag", 'b', omitable = true);
  parameter -> SetDefaultValue (true);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("axes-flag", 'b', omitable = true);
  parameter -> SetDefaultValue (true);
  parameter -> SetGuidance ("Set \"false\" to suppress axes.");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("check-overlap-flag", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  parameter -> SetGuidance ("Set \"false\" to suppress overlap check.");
  fpCommand->SetParameter(parameter);
}

G4VisCommandSceneAddLogicalVolume::~G4VisCommandSceneAddLogicalVolume () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddLogicalVolume::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddLogicalVolume::SetNewValue (G4UIcommand*,
						     G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String name;
  G4int requestedDepthOfDescent;
  G4String booleansString, voxelsString, readoutString, axesString;
  G4String overlapString;
  std::istringstream is (newValue);
  is >> name >> requestedDepthOfDescent
     >>  booleansString >> voxelsString >> readoutString >> axesString
     >> overlapString;
  G4bool booleans = G4UIcommand::ConvertToBool(booleansString);
  G4bool voxels = G4UIcommand::ConvertToBool(voxelsString);
  G4bool readout = G4UIcommand::ConvertToBool(readoutString);
  G4bool axes = G4UIcommand::ConvertToBool(axesString);
  G4bool checkOverlaps = G4UIcommand::ConvertToBool(overlapString);

  G4LogicalVolumeStore *pLVStore = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* pLV = nullptr;
  pLV = pLVStore->GetVolume(name);
  if (pLV == nullptr) return;  // Volume not found; warning message thrown

  const std::vector<G4Scene::Model>& rdModelList =
    pScene -> GetRunDurationModelList();
  std::vector<G4Scene::Model>::const_iterator i;
  for (i = rdModelList.begin(); i != rdModelList.end(); ++i) {
    if (i->fpModel->GetGlobalDescription().find("Volume") != std::string::npos) break;
  }
  if (i != rdModelList.end()) {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "There is already a volume, \""
             << i->fpModel->GetGlobalDescription()
             << "\",\n in the run-duration model list of scene \""
             << pScene -> GetName()
             << "\".\n Your logical volume must be the only volume in the scene."
	     << "\n Create a new scene and try again:"
	     << "\n  /vis/specify " << name
	     << "\n or"
	     << "\n  /vis/scene/create"
	     << "\n  /vis/scene/add/logicalVolume " << name
	     << "\n  /vis/sceneHandler/attach"
	     << "\n (and also, if necessary, /vis/viewer/flush)"
             << G4endl;
    }
    return;
  }

  G4LogicalVolumeModel* model = new G4LogicalVolumeModel
    (pLV, requestedDepthOfDescent, booleans, voxels, readout, checkOverlaps);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);

  if (successful) {

    G4bool axesSuccessful = false;
    if (axes) {
      const G4double radius = model->GetExtent().GetExtentRadius();
      const G4double axisLengthMax = radius / 2.;
      const G4double intLog10Length = std::floor(std::log10(axisLengthMax));
      G4double axisLength = std::pow(10,intLog10Length);
      if (5.*axisLength < axisLengthMax) axisLength *= 5.;
      else if (2.*axisLength < axisLengthMax) axisLength *= 2.;
      const G4double axisWidth = axisLength / 20.;
      G4VModel* axesModel = new G4AxesModel(0.,0.,0.,axisLength,axisWidth);
      axesSuccessful = pScene -> AddRunDurationModel (axesModel, warn);
    }

//    if (verbosity >= G4VisManager::warnings) {
//      const std::map<G4String,G4AttDef>* attDefs = model->GetAttDefs();
//      std::vector<G4AttValue>* attValues = model->CreateCurrentAttValues();
//      G4warn << G4AttCheck(attValues, attDefs);
//      delete attValues;
//    }

    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Logical volume \"" << pLV -> GetName ()
	     << "\" with requested depth of descent "
	     << requestedDepthOfDescent
	     << ",\n  with";
      if (!booleans) G4cout << "out";
      G4cout << " boolean components, with";
      if (!voxels) G4cout << "out";
      G4cout << " voxels,\n  with";
      if (!readout) G4cout << "out";
      G4cout << " readout geometry and with";
      if (!checkOverlaps) G4cout << "out";
      G4cout << " overlap checking"
	     << "\n  has been added to scene \"" << currentSceneName << "\".";
      if (axes) {
        if (axesSuccessful) {
          G4cout <<
          "\n  Axes have also been added at the origin of local cooordinates.";
        } else {
          G4cout <<
          "\n  Axes have not been added for some reason possibly stated above.";
        }
      }
      G4cout << G4endl;
    }
  }
  else {
    G4VisCommandsSceneAddUnsuccessful(verbosity);
    return;
  }

  CheckSceneAndNotifyHandlers (pScene);
}


////////////// /vis/scene/add/logo //////////////////////////////////

G4VisCommandSceneAddLogo::G4VisCommandSceneAddLogo () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/logo", this);
  fpCommand -> SetGuidance ("Adds a G4 logo to the current scene.");
  fpCommand -> SetGuidance
  ("If \"unit\" is \"auto\", height is roughly one tenth of scene extent.");
  fpCommand -> SetGuidance
  ("\"direction\" is that of outward-facing normal to front face of logo."
   "\nIf \"direction\" is \"auto\", logo faces the user in the current viewer.");
  fpCommand -> SetGuidance
  ("\nIf \"placement\" is \"auto\", logo is placed at bottom right of screen"
   "\n  when viewed from logo direction.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("height", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("auto");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("direction", 's', omitable = true);
  parameter->SetGuidance ("auto|[-]x|[-]y|[-]z");
  parameter->SetDefaultValue ("auto");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("red", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("placement", 's', omitable = true);
  parameter -> SetParameterCandidates("auto manual");
  parameter->SetDefaultValue  ("auto");
  fpCommand->SetParameter     (parameter);
  parameter =  new G4UIparameter ("xmid", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("ymid", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("zmid", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSceneAddLogo::~G4VisCommandSceneAddLogo () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddLogo::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddLogo::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn = verbosity >= G4VisManager::warnings;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  } else {
    if (pScene->GetExtent().GetExtentRadius() <= 0.) {
      if (verbosity >= G4VisManager::errors) {
        G4warn
  << "ERROR: Scene has no extent. Add volumes or use \"/vis/scene/add/extent\"."
        << G4endl;
      }
      return;
    }
  }

  G4VViewer* pViewer = fpVisManager->GetCurrentViewer();
  if (!pViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
	"ERROR: G4VisCommandSceneAddLogo::SetNewValue: no viewer."
	"\n  Auto direction needs a viewer."
	     << G4endl;
    }
    return;
  }

  G4double userHeight, red, green, blue, xmid, ymid, zmid;
  G4String userHeightUnit, direction, placement, positionUnit;
  std::istringstream is (newValue);
  is >> userHeight >> userHeightUnit >> direction
     >> red >> green >> blue
     >> placement
     >> xmid >> ymid >> zmid >> positionUnit;

  G4double height = userHeight;
  const G4VisExtent& sceneExtent = pScene->GetExtent();  // Existing extent.
  if (userHeightUnit == "auto") {
    height *= 0.2 * sceneExtent.GetExtentRadius();
  } else {
    height *= G4UIcommand::ValueOf(userHeightUnit);
  }

  G4double unit = G4UIcommand::ValueOf(positionUnit);
  xmid *= unit; ymid *= unit; zmid *= unit;

  Direction logoDirection = X;  // Initialise to keep some compilers happy.
  if (direction == "auto") {
    // Take cue from viewer
    const G4Vector3D& vp =
      pViewer->GetViewParameters().GetViewpointDirection();
    if (vp.x() > vp.y() && vp.x() > vp.z()) logoDirection = X;
    else if (vp.x() < vp.y() && vp.x() < vp.z()) logoDirection = minusX;
    else if (vp.y() > vp.x() && vp.y() > vp.z()) logoDirection = Y;
    else if (vp.y() < vp.x() && vp.y() < vp.z()) logoDirection = minusY;
    else if (vp.z() > vp.x() && vp.z() > vp.y()) logoDirection = Z;
    else if (vp.z() < vp.x() && vp.z() < vp.y()) logoDirection = minusZ;
  }
  else if (direction[0] == 'x') logoDirection = X;
  else if (direction[0] == 'y') logoDirection = Y;
  else if (direction[0] == 'z') logoDirection = Z;
  else if (direction[0] == '-') {
    if (direction[1] == 'x') logoDirection = minusX;
    else if (direction[1] == 'y') logoDirection = minusY;
    else if (direction[1] == 'z') logoDirection = minusZ;
  } else {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: Unrecogniseed direction: \""
	     << direction << "\"." << G4endl;
      return;
    }
  }

  G4bool autoPlacing = false; if (placement == "auto") autoPlacing = true;
  // Parameters read and interpreted.

  // Current scene extent
  const G4double xmin = sceneExtent.GetXmin();
  const G4double xmax = sceneExtent.GetXmax();
  const G4double ymin = sceneExtent.GetYmin();
  const G4double ymax = sceneExtent.GetYmax();
  const G4double zmin = sceneExtent.GetZmin();
  const G4double zmax = sceneExtent.GetZmax();

  // Test existing extent and issue warnings...
  G4bool worried = false;
  if (sceneExtent.GetExtentRadius() == 0) {
    worried = true;
    if (verbosity >= G4VisManager::warnings) {
      G4warn <<
	"WARNING: Existing scene does not yet have any extent."
	"\n  Maybe you have not yet added any geometrical object."
	     << G4endl;
    }
  }

  // Useful constants, etc...
  const G4double halfHeight(height / 2.);
  const G4double comfort(0.01);  // 0.15 seems too big.  0.05 might be better.
  const G4double freeHeightFraction (1. + 2. * comfort);

  // Test existing scene for room...
  G4bool room = true;
  switch (logoDirection) {
  case X:
  case minusX:
    if (freeHeightFraction * (xmax - xmin) < height) room = false;
    break;
  case Y:
  case minusY:
    if (freeHeightFraction * (ymax - ymin) < height) room = false;
    break;
  case Z:
  case minusZ:
    if (freeHeightFraction * (zmax - zmin) < height) room = false;
    break;
  }
  if (!room) {
    worried = true;
    if (verbosity >= G4VisManager::warnings) {
      G4warn <<
	"WARNING: Not enough room in existing scene.  Maybe logo is too large."
	     << G4endl;
    }
  }
  if (worried) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn <<
	"WARNING: The logo you have asked for is bigger than the existing"
	"\n  scene.  Maybe you have added it too soon.  It is recommended that"
	"\n  you add the logo last so that it can be correctly auto-positioned"
	"\n  so as not to be obscured by any existing object and so that the"
	"\n  view parameters can be correctly recalculated."
	     << G4endl;
    }
  }

  G4double sxmid(xmid), symid(ymid), szmid(zmid);
  if (autoPlacing) {
    // Aim to place at bottom right of screen when viewed from logoDirection.
    // Give some comfort zone.
    const G4double xComfort = comfort * (xmax - xmin);
    const G4double yComfort = comfort * (ymax - ymin);
    const G4double zComfort = comfort * (zmax - zmin);
    switch (logoDirection) {
    case X:  // y-axis up, z-axis to left?
      sxmid = xmax + halfHeight + xComfort;
      symid = ymin - yComfort;
      szmid = zmin - zComfort;
      break;
    case minusX:  // y-axis up, z-axis to right?
      sxmid = xmin - halfHeight - xComfort;
      symid = ymin - yComfort;
      szmid = zmax + zComfort;
      break;
    case Y:  // z-axis up, x-axis to left?
      sxmid = xmin - xComfort;
      symid = ymax + halfHeight + yComfort;
      szmid = zmin - zComfort;
      break;
    case minusY:  // z-axis up, x-axis to right?
      sxmid = xmax + xComfort;
      symid = ymin - halfHeight - yComfort;
      szmid = zmin - zComfort;
      break;
    case Z:  // y-axis up, x-axis to right?
      sxmid = xmax + xComfort;
      symid = ymin - yComfort;
      szmid = zmax + halfHeight + zComfort;
      break;
    case minusZ:  // y-axis up, x-axis to left?
      sxmid = xmin - xComfort;
      symid = ymin - yComfort;
      szmid = zmin - halfHeight - zComfort;
      break;
    }
  }

  G4Transform3D transform;
  switch (logoDirection) {
  case X:  // y-axis up, z-axis to left?
    transform = G4RotateY3D(halfpi);
    break;
  case minusX:  // y-axis up, z-axis to right?
    transform = G4RotateY3D(-halfpi);
    break;
  case Y:  // z-axis up, x-axis to left?
    transform = G4RotateX3D(-halfpi) * G4RotateZ3D(pi);
    break;
  case minusY:  // z-axis up, x-axis to right?
    transform = G4RotateX3D(halfpi);
    break;
  case Z:  // y-axis up, x-axis to right?
    // No transformation required.
    break;
  case minusZ:  // y-axis up, x-axis to left?
    transform = G4RotateY3D(pi);
    break;
  }
  transform = G4Translate3D(sxmid,symid,szmid) * transform;

  G4VisAttributes visAtts(G4Colour(red, green, blue));
  visAtts.SetForceSolid(true);         // Always solid.

  G4Logo* logo = new G4Logo(height,visAtts,transform);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddLogo::G4Logo>(logo);
  model->SetType("G4Logo");
  model->SetGlobalTag("G4Logo");
  model->SetGlobalDescription("G4Logo: " + newValue);
  G4double& h = height;
  G4double h2 = h/2.;
  G4VisExtent extent(-h,h,-h2,h2,-h2,h2);
  model->SetExtent(extent.Transform(transform));
  // This extent gets "added" to existing scene extent in
  // AddRunDurationModel below.
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "G4 Logo of height " << userHeight << ' ' << userHeightUnit
	     << ", " << direction << "-direction, added to scene \""
	     << currentSceneName << "\"";
      if (verbosity >= G4VisManager::parameters) {
	G4cout << "\n  with extent " << extent
	       << "\n  at " << transform.getRotation()
	       << "  " << transform.getTranslation();
      }
      G4cout << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

G4VisCommandSceneAddLogo::G4Logo::G4Logo
(G4double height, const G4VisAttributes& visAtts, const G4Transform3D& transform)
{
  const G4double& h =  height;
  const G4double h2  = 0.5 * h;   // Half height.
  const G4double ri  = 0.25 * h;  // Inner radius.
  const G4double ro  = 0.5 * h;   // Outer radius.
  const G4double ro2 = 0.5 * ro;  // Half outer radius.
  const G4double w   = ro - ri;   // Width.
  const G4double w2  = 0.5 * w;   // Half width.
  const G4double d2  = 0.2 * h;   // Half depth.
  const G4double f1  = 0.05 * h;  // left edge of stem of "4".
  const G4double f2  = -0.3 * h;  // bottom edge of cross of "4".
  const G4double e = 1.e-4 * h;   // epsilon.
  const G4double xt = f1, yt = h2;      // Top of slope.
  const G4double xb = -h2, yb = f2 + w; // Bottom of slope.
  const G4double dx = xt - xb, dy = yt - yb;
  const G4double angle = std::atan2(dy,dx);
  G4RotationMatrix rm;
  rm.rotateZ(angle*rad);
  const G4double d = std::sqrt(dx * dx + dy * dy);
  const G4double ss = h;  // Half height of square subtractor
  const G4double y8 = ss; // Choose y of subtractor for outer slope.
  const G4double x8 = ((-ss * d - dx * (yt - y8)) / dy) + xt;
  G4double y9 = ss; // Choose y of subtractor for inner slope.
  G4double x9 = ((-(ss - w) * d - dx * (yt - y8)) / dy) + xt;
  // But to get inner, we make a triangle translated by...
  const G4double xtr = ss - f1, ytr = -ss - f2 -w;
  x9 += xtr; y9 += ytr;

  // The idea here is to create a polyhedron for the G and the 4.  To do
  // this we use Geant4 geometry solids and make boolean operations.
  // Note that we do not need to keep the solids. We use local objects,
  // which, of course, are deleted on leaving this function. This
  // is contrary to the usual requirement for solids that are part of the
  // detector for which solids MUST be created on the heap (with "new").
  // Finally we invoke CreatePolyhedron, which creates a polyhedron on the heap
  // and returns a pointer.  It is the user's responsibility to delete,
  // which is done in the destructor of this class. Thus the polyhedra,
  // created here, remain on the heap for the lifetime of the job.

  // G...
  G4Tubs tG("tG",ri,ro,d2,0.15*pi,1.85*pi);
  G4Box bG("bG",w2,ro2,d2);
  G4UnionSolid logoG("logoG",&tG,&bG,G4Translate3D(ri+w2,-ro2,0.));
  fpG = logoG.CreatePolyhedron();
  fpG->SetVisAttributes(visAtts);
  fpG->Transform(G4Translate3D(-0.55*h,0.,0.));
  fpG->Transform(transform);

  // 4...
  G4Box b1("b1",h2,h2,d2);
  G4Box bS("bS",ss,ss,d2+e);  // Subtractor.
  G4Box bS2("bS2",ss,ss,d2+2.*e);  // 2nd Subtractor.
  G4SubtractionSolid s1("s1",&b1,&bS,G4Translate3D(f1-ss,f2-ss,0.));
  G4SubtractionSolid s2("s2",&s1,&bS,G4Translate3D(f1+ss+w,f2-ss,0.));
  G4SubtractionSolid s3("s3",&s2,&bS,G4Translate3D(f1+ss+w,f2+ss+w,0.));
  G4SubtractionSolid s4
    ("s4",&s3,&bS,G4Transform3D(rm,G4ThreeVector(x8,y8,0.)));
  G4SubtractionSolid s5    // Triangular hole.
    ("s5",&bS,&bS2,G4Transform3D(rm,G4ThreeVector(x9,y9,0.)));
  G4SubtractionSolid logo4("logo4",&s4,&s5,G4Translate3D(-xtr,-ytr,0.));
  fp4 = logo4.CreatePolyhedron();
  /* Experiment with creating own polyhedron...
  int nNodes = 4;
  int nFaces = 4;
  double xyz[][3] = {{0,0,0},{1*m,0,0},{0,1*m,0},{0,0,1*m}};
  int faces[][4] = {{1,3,2,0},{1,2,4,0},{1,4,3,0},{2,3,4,0}};
  fp4 = new G4Polyhedron();
  fp4->createPolyhedron(nNodes,nFaces,xyz,faces);
  */
  fp4->SetVisAttributes(visAtts);
  fp4->Transform(G4Translate3D(0.55*h,0.,0.));
  fp4->Transform(transform);
}

G4VisCommandSceneAddLogo::G4Logo::~G4Logo() {
  delete fpG;
  delete fp4;
}

void G4VisCommandSceneAddLogo::G4Logo::operator()
  (G4VGraphicsScene& sceneHandler, const G4ModelingParameters*) {
  sceneHandler.BeginPrimitives();
  sceneHandler.AddPrimitive(*fpG);
  sceneHandler.AddPrimitive(*fp4);
  sceneHandler.EndPrimitives();
}

////////////// /vis/scene/add/logo2D ///////////////////////////////////////

G4VisCommandSceneAddLogo2D::G4VisCommandSceneAddLogo2D () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/logo2D", this);
  fpCommand -> SetGuidance ("Adds 2D logo to current scene.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("size", 'i', omitable = true);
  parameter -> SetGuidance ("Screen size of text in pixels.");
  parameter -> SetDefaultValue (48);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("x-position", 'd', omitable = true);
  parameter -> SetGuidance ("x screen position in range -1 < x < 1.");
  parameter -> SetDefaultValue (-0.9);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("y-position", 'd', omitable = true);
  parameter -> SetGuidance ("y screen position in range -1 < y < 1.");
  parameter -> SetDefaultValue (-0.9);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("layout", 's', omitable = true);
  parameter -> SetGuidance ("Layout, i.e., adjustment: left|centre|right.");
  parameter -> SetDefaultValue ("left");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddLogo2D::~G4VisCommandSceneAddLogo2D () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddLogo2D::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddLogo2D::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int size;
  G4double x, y;
  G4String layoutString;
  std::istringstream is(newValue);
  is >> size >> x >> y >> layoutString;
  G4Text::Layout layout = G4Text::right;
  if (layoutString[0] == 'l') layout = G4Text::left;
  else if (layoutString[0] == 'c') layout = G4Text::centre;
  else if (layoutString[0] == 'r') layout = G4Text::right;

  Logo2D* logo2D = new Logo2D(fpVisManager, size, x, y, layout);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddLogo2D::Logo2D>(logo2D);
  model->SetType("G4Logo2D");
  model->SetGlobalTag("G4Logo2D");
  model->SetGlobalDescription("G4Logo2D: " + newValue);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "2D logo has been added to scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

void G4VisCommandSceneAddLogo2D::Logo2D::operator()
  (G4VGraphicsScene& sceneHandler, const G4ModelingParameters*)
{
  G4Text text("Geant4", G4Point3D(fX, fY, 0.));
  text.SetScreenSize(fSize);
  text.SetLayout(fLayout);
  G4VisAttributes textAtts(G4Colour::Brown());
  text.SetVisAttributes(textAtts);
  sceneHandler.BeginPrimitives2D();
  sceneHandler.AddPrimitive(text);
  sceneHandler.EndPrimitives2D();
}

////////////// /vis/scene/add/magneticField ///////////////////////////////////////

G4VisCommandSceneAddMagneticField::G4VisCommandSceneAddMagneticField () {
  fpCommand = new G4UIcommand ("/vis/scene/add/magneticField", this);
  fpCommand -> SetGuidance
  ("Adds magnetic field representation to current scene.");
  const G4UIcommandTree* tree = G4UImanager::GetUIpointer()->GetTree();
  const G4UIcommand* addElecFieldCmd = tree->FindPath("/vis/scene/add/electricField");
  // Pick up additional guidance from /vis/scene/add/electricField
  CopyGuidanceFrom(addElecFieldCmd,fpCommand,1);
  // Pick up parameters from /vis/scene/add/electricField
  CopyParametersFrom(addElecFieldCmd,fpCommand);
}

G4VisCommandSceneAddMagneticField::~G4VisCommandSceneAddMagneticField () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddMagneticField::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddMagneticField::SetNewValue
(G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<  "ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int nDataPointsPerHalfScene;
  G4String representation;
  std::istringstream iss(newValue);
  iss >> nDataPointsPerHalfScene >> representation;
  G4ElectricFieldModel::Representation
  modelRepresentation = G4ElectricFieldModel::fullArrow;
  if (representation == "lightArrow") {
    modelRepresentation = G4ElectricFieldModel::lightArrow;
  }
  G4VModel* model;
  model = new G4MagneticFieldModel
  (nDataPointsPerHalfScene,modelRepresentation,
   fCurrentArrow3DLineSegmentsPerCircle,
   fCurrentExtentForField,
   fCurrrentPVFindingsForField);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout
      << "Magnetic field, if any, will be drawn in scene \""
      << currentSceneName
      << "\"\n  with "
      << nDataPointsPerHalfScene
      << " data points per half extent and with representation \""
      << representation
      << '\"'
      << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/add/psHits ///////////////////////////////////////

G4VisCommandSceneAddPSHits::G4VisCommandSceneAddPSHits () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/add/psHits", this);
  fpCommand -> SetGuidance
    ("Adds Primitive Scorer Hits (PSHits) to current scene.");
  fpCommand -> SetGuidance
    ("PSHits are drawn at end of run when the scene in which"
     "\nthey are added is current.");
  fpCommand -> SetGuidance
    ("Optional parameter specifies name of scoring map.  By default all"
     "\nscoring maps registered with the G4ScoringManager are drawn.");
  fpCommand -> SetParameterName ("mapname", omitable = true);
  fpCommand -> SetDefaultValue ("all");
}

G4VisCommandSceneAddPSHits::~G4VisCommandSceneAddPSHits () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddPSHits::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddPSHits::SetNewValue
(G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4VModel* model = new G4PSHitsModel(newValue);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddEndOfRunModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      if (newValue == "all") {
	G4cout << "All Primitive Scorer hits";
      } else {
	G4cout << "Hits of Primitive Scorer \"" << newValue << '"';
      }
      G4cout << " will be drawn at end of run in scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/add/scale //////////////////////////////////

G4VisCommandSceneAddScale::G4VisCommandSceneAddScale () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/scale", this);
  fpCommand -> SetGuidance
  ("Adds an annotated scale line to the current scene.");
  fpCommand -> SetGuidance
  ("If \"unit\" is \"auto\", length is roughly one tenth of the scene extent.");
  fpCommand -> SetGuidance
  ("If \"direction\" is \"auto\", scale is roughly in the plane of the current view.");
  fpCommand -> SetGuidance
  ("If \"placement\" is \"auto\", scale is placed at bottom left of current view."
   "\n  Otherwise placed at (xmid,ymid,zmid).");
  fpCommand -> SetGuidance
  ("An annotated line in the specified direction with tick marks at the"
   "\nend.  If autoPlacing is true it is required to be centred at the"
   "\nfront, right, bottom corner of the world space, comfortably outside"
   "\nthe existing bounding box/sphere so that existing objects do not"
   "\nobscure it.  Otherwise it is required to be drawn with mid-point at"
   "\n(xmid, ymid, zmid)."
   "\n"
   "\nThe auto placing algorithm is (approx):"
   "\n  x = xmin + (1 + comfort) * (xmax - xmin);"
   "\n  y = ymin - comfort * (ymax - ymin);"
   "\n  z = zmin + (1 + comfort) * (zmax - zmin);"
   "\n  if direction == x then (x - length,y,z) to (x,y,z);"
   "\n  if direction == y then (x,y,z) to (x,y + length,z);"
   "\n  if direction == z then (x,y,z - length) to (x,y,z);"
   );
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("length", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("auto");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("direction", 's', omitable = true);
  parameter->SetGuidance ("auto|x|y|z");
  parameter->SetDefaultValue ("auto");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("red", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("placement", 's', omitable = true);
  parameter -> SetParameterCandidates("auto manual");
  parameter->SetDefaultValue  ("auto");
  fpCommand->SetParameter     (parameter);
  parameter =  new G4UIparameter ("xmid", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("ymid", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("zmid", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSceneAddScale::~G4VisCommandSceneAddScale () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddScale::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddScale::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn = verbosity >= G4VisManager::warnings;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  } else {
    if (pScene->GetExtent().GetExtentRadius() <= 0.) {
      if (verbosity >= G4VisManager::errors) {
        G4warn
  << "ERROR: Scene has no extent. Add volumes or use \"/vis/scene/add/extent\"."
        << G4endl;
      }
      return;
    }
  }

  G4double userLength, red, green, blue, xmid, ymid, zmid;
  G4String userLengthUnit, direction, placement, positionUnit;
  std::istringstream is (newValue);
  is >> userLength >> userLengthUnit >> direction
     >> red >> green >> blue
     >> placement
     >> xmid >> ymid >> zmid >> positionUnit;

  G4double length = userLength;
  const G4VisExtent& sceneExtent = pScene->GetExtent();  // Existing extent.
  if (userLengthUnit == "auto") {
    const G4double lengthMax = 0.5 * sceneExtent.GetExtentRadius();
    const G4double intLog10Length = std::floor(std::log10(lengthMax));
    length = std::pow(10,intLog10Length);
    if (5.*length < lengthMax) length *= 5.;
    else if (2.*length < lengthMax) length *= 2.;
  } else {
    length *= G4UIcommand::ValueOf(userLengthUnit);
  }
  G4String annotation = G4BestUnit(length,"Length");

  G4double unit = G4UIcommand::ValueOf(positionUnit);
  xmid *= unit; ymid *= unit; zmid *= unit;

  Scale::Direction scaleDirection (Scale::x);
  if (direction[0] == 'y') scaleDirection = Scale::y;
  if (direction[0] == 'z') scaleDirection = Scale::z;

  G4VViewer* pViewer = fpVisManager->GetCurrentViewer();
  if (!pViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
	"ERROR: G4VisCommandSceneAddScale::SetNewValue: no viewer."
	"\n  Auto direction needs a viewer."
	     << G4endl;
    }
    return;
  }

  const G4Vector3D& vp =
    pViewer->GetViewParameters().GetViewpointDirection();
  const G4Vector3D& up =
    pViewer->GetViewParameters().GetUpVector();

  if (direction == "auto") {  // Takes cue from viewer.
    if (std::abs(vp.x()) > std::abs(vp.y()) &&
	std::abs(vp.x()) > std::abs(vp.z())) {  // x viewpoint
      if (std::abs(up.y()) > std::abs(up.z())) scaleDirection = Scale::z;
	  else scaleDirection = Scale::y;
    }
    else if (std::abs(vp.y()) > std::abs(vp.x()) &&
	     std::abs(vp.y()) > std::abs(vp.z())) {  // y viewpoint
      if (std::abs(up.x()) > std::abs(up.z())) scaleDirection = Scale::z;
	  else scaleDirection = Scale::x;
    }
    else if (std::abs(vp.z()) > std::abs(vp.x()) &&
	     std::abs(vp.z()) > std::abs(vp.y())) {  // z viewpoint
      if (std::abs(up.y()) > std::abs(up.x())) scaleDirection = Scale::x;
	  else scaleDirection = Scale::y;
    }
  }

  G4bool autoPlacing = false; if (placement == "auto") autoPlacing = true;
  // Parameters read and interpreted.

  // Useful constants, etc...
  const G4double halfLength(length / 2.);
  const G4double comfort(0.01);  // 0.15 seems too big.  0.05 might be better.
  const G4double freeLengthFraction (1. + 2. * comfort);

  const G4double xmin = sceneExtent.GetXmin();
  const G4double xmax = sceneExtent.GetXmax();
  const G4double ymin = sceneExtent.GetYmin();
  const G4double ymax = sceneExtent.GetYmax();
  const G4double zmin = sceneExtent.GetZmin();
  const G4double zmax = sceneExtent.GetZmax();

  // Test existing extent and issue warnings...
  G4bool worried = false;
  if (sceneExtent.GetExtentRadius() == 0) {
    worried = true;
    if (verbosity >= G4VisManager::warnings) {
      G4warn <<
	"WARNING: Existing scene does not yet have any extent."
	"\n  Maybe you have not yet added any geometrical object."
	     << G4endl;
    }
  }

  // Test existing scene for room...
  G4bool room  = true;
  switch (scaleDirection) {
  case Scale::x:
    if (freeLengthFraction * (xmax - xmin) < length) room = false;
    break;
  case Scale::y:
    if (freeLengthFraction * (ymax - ymin) < length) room = false;
    break;
  case Scale::z:
    if (freeLengthFraction * (zmax - zmin) < length) room = false;
    break;
  }
  if (!room) {
    worried = true;
    if (verbosity >= G4VisManager::warnings) {
      G4warn <<
	"WARNING: Not enough room in existing scene.  Maybe scale is too long."
	     << G4endl;
    }
  }
  if (worried) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn <<
	"WARNING: The scale you have asked for is bigger than the existing"
	"\n  scene.  Maybe you have added it too soon.  It is recommended that"
	"\n  you add the scale last so that it can be correctly auto-positioned"
	"\n  so as not to be obscured by any existing object and so that the"
	"\n  view parameters can be correctly recalculated."
	     << G4endl;
    }
  }

  // Now figure out the extent...
  //
  // This creates a representation of annotated line in the specified
  // direction with tick marks at the end.  If autoPlacing is true it
  // is required to be centred at the front, right, bottom corner of
  // the world space, comfortably outside the existing bounding
  // box/sphere so that existing objects do not obscure it.  Otherwise
  // it is required to be drawn with mid-point at (xmid, ymid, zmid).
  //
  // The auto placing algorithm might be:
  //   x = xmin + (1 + comfort) * (xmax - xmin)
  //   y = ymin - comfort * (ymax - ymin)
  //   z = zmin + (1 + comfort) * (zmax - zmin)
  //   if direction == x then (x - length,y,z) to (x,y,z)
  //   if direction == y then (x,y,z) to (x,y + length,z)
  //   if direction == z then (x,y,z - length) to (x,y,z)
  //
  // Implement this in two parts.  Here, use the scale's extent to
  // "expand" the scene's extent.  Then rendering - in
  // G4VSceneHandler::AddPrimitive(const G4Scale&) - simply has to
  // ensure it's within the new extent.
  //

  G4double sxmid(xmid), symid(ymid), szmid(zmid);
  if (autoPlacing) {
    // Aim to place at bottom right of screen in current view.
    // Give some comfort zone.
    const G4double xComfort = comfort * (xmax - xmin);
    const G4double yComfort = comfort * (ymax - ymin);
    const G4double zComfort = comfort * (zmax - zmin);
    switch (scaleDirection) {
    case Scale::x:
      if (vp.z() > 0.) {
	sxmid = xmax + xComfort;
	symid = ymin - yComfort;
	szmid = zmin - zComfort;
      } else {
	sxmid = xmin - xComfort;
	symid = ymin - yComfort;
	szmid = zmax + zComfort;
      }
      break;
    case Scale::y:
      if (vp.x() > 0.) {
	sxmid = xmin - xComfort;
	symid = ymax + yComfort;
	szmid = zmin - zComfort;
      } else {
	sxmid = xmax + xComfort;
	symid = ymin - yComfort;
	szmid = zmin - zComfort;
      }
      break;
    case Scale::z:
      if (vp.x() > 0.) {
	sxmid = xmax + xComfort;
	symid = ymin - yComfort;
	szmid = zmax + zComfort;
      } else {
	sxmid = xmin - xComfort;
	symid = ymin - yComfort;
	szmid = zmax + zComfort;
      }
      break;
    }
  }

  G4Transform3D transform;
  const G4double h = halfLength;
  const G4double t = h/5.;
  G4VisExtent scaleExtent(-h,h,-t,t,-t,t);
  switch (scaleDirection) {
  case Scale::x:
    break;
  case Scale::y:
    transform = G4RotateZ3D(halfpi);
    break;
  case Scale::z:
    transform = G4RotateY3D(halfpi);
    break;
  }
  transform = G4Translate3D(sxmid,symid,szmid) * transform;
  scaleExtent = scaleExtent.Transform(transform);

  G4Colour colour(red, green, blue);
  if (direction == "auto") {
    switch (scaleDirection) {
      case Scale::x:
	colour = G4Colour::Red();
	break;
      case Scale::y:
	colour = G4Colour::Green();
	break;
      case Scale::z:
	colour = G4Colour::Blue();
	break;
    }
  }
  G4VisAttributes visAttr(colour);

  Scale* scale = new Scale
  (visAttr, length, transform,
   annotation, fCurrentTextSize, colour);
  G4VModel* model = new G4CallbackModel<Scale>(scale);
  model->SetType("Scale");
  model->SetGlobalTag("Scale");
  model->SetGlobalDescription("Scale: " + newValue);
  model->SetExtent(scaleExtent);

  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Scale of " << annotation
	     << " added to scene \"" << currentSceneName << "\".";
      if (verbosity >= G4VisManager::parameters) {
	G4cout << "\n  with extent " << scaleExtent
	       << "\n  at " << transform.getRotation()
	       << "  " << transform.getTranslation();
      }
      G4cout << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

G4VisCommandSceneAddScale::Scale::Scale
 (const G4VisAttributes& visAtts,
  G4double length, const G4Transform3D& transform,
  const G4String& annotation, G4double annotationSize,
  const G4Colour& annotationColour):
fVisAtts(visAtts)
{
  // Useful constants...
  const G4double halfLength(length / 2.);
  const G4double tickLength(length / 20.);

  // Create (empty) polylines having the same vis attributes...
  // (OK to pass address since fVisAtts is long lived.)
  fScaleLine.SetVisAttributes(&fVisAtts);
  fTick11.SetVisAttributes(&fVisAtts);
  fTick12.SetVisAttributes(&fVisAtts);
  fTick21.SetVisAttributes(&fVisAtts);
  fTick22.SetVisAttributes(&fVisAtts);

  // Add points to the polylines to represent a scale parallel to the
  // x-axis centred on the origin...
  G4Point3D r1(G4Point3D(-halfLength, 0., 0.));
  G4Point3D r2(G4Point3D( halfLength, 0., 0.));
  fScaleLine.push_back(r1);
  fScaleLine.push_back(r2);
  G4Point3D ticky(0., tickLength, 0.);
  G4Point3D tickz(0., 0., tickLength);
  fTick11.push_back(r1 + ticky);
  fTick11.push_back(r1 - ticky);
  fTick12.push_back(r1 + tickz);
  fTick12.push_back(r1 - tickz);
  fTick21.push_back(r2 + ticky);
  fTick21.push_back(r2 - ticky);
  fTick22.push_back(r2 + tickz);
  fTick22.push_back(r2 - tickz);
  // ...and transform to chosen position and orientation
  fScaleLine.transform(transform);
  fTick11.transform(transform);
  fTick12.transform(transform);
  fTick21.transform(transform);
  fTick22.transform(transform);
  // Similarly for annotation
  G4Point3D textPosition(0., tickLength, 0.);
  textPosition.transform(transform);
  fText = G4Text(annotation,textPosition);
  fText.SetVisAttributes(annotationColour);
  fText.SetScreenSize(annotationSize);
}

void G4VisCommandSceneAddScale::Scale::operator()
(G4VGraphicsScene& sceneHandler,const G4ModelingParameters*)
{
  // Draw...
  sceneHandler.BeginPrimitives();
  sceneHandler.AddPrimitive(fScaleLine);
  sceneHandler.AddPrimitive(fTick11);
  sceneHandler.AddPrimitive(fTick12);
  sceneHandler.AddPrimitive(fTick21);
  sceneHandler.AddPrimitive(fTick22);
  sceneHandler.AddPrimitive(fText);
  sceneHandler.EndPrimitives();
}

////////////// /vis/scene/add/text //////////////////////////////////

G4VisCommandSceneAddText::G4VisCommandSceneAddText () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/text", this);
  fpCommand -> SetGuidance ("Adds text to current scene.");
  fpCommand -> SetGuidance
    ("Use \"/vis/set/textColour\" to set colour.");
  fpCommand -> SetGuidance
    ("Use \"/vis/set/textLayout\" to set layout:");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("x", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("y", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("z", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter     (parameter);
  parameter =  new G4UIparameter ("font_size", 'd', omitable = true);
  parameter->SetDefaultValue (12);
  parameter->SetGuidance ("pixels");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("x_offset", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  parameter->SetGuidance ("pixels");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("y_offset", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  parameter->SetGuidance ("pixels");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("text", 's', omitable = true);
  parameter->SetGuidance ("The rest of the line is text.");
  parameter->SetDefaultValue ("Hello G4");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSceneAddText::~G4VisCommandSceneAddText () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddText::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddText::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn = verbosity >= G4VisManager::warnings;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4Tokenizer next(newValue);
  G4double x = StoD(next());
  G4double y = StoD(next());
  G4double z = StoD(next());
  G4String unitString = next();
  G4double font_size = StoD(next());
  G4double x_offset = StoD(next());
  G4double y_offset = StoD(next());
  G4String text = next("\n");

  G4double unit = G4UIcommand::ValueOf(unitString);
  x *= unit; y *= unit; z *= unit;

  G4Text g4text(text, G4Point3D(x,y,z));
  G4VisAttributes visAtts(fCurrentTextColour);
  g4text.SetVisAttributes(visAtts);
  g4text.SetLayout(fCurrentTextLayout);
  g4text.SetScreenSize(font_size);
  g4text.SetOffset(x_offset,y_offset);
  G4VModel* model = new G4TextModel(g4text);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Text \"" << text
	     << "\" has been added to scene \"" << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}


////////////// /vis/scene/add/text2D //////////////////////////////////

G4VisCommandSceneAddText2D::G4VisCommandSceneAddText2D () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/text2D", this);
  fpCommand -> SetGuidance ("Adds 2D text to current scene.");
  fpCommand -> SetGuidance
    ("Use \"/vis/set/textColour\" to set colour.");
  fpCommand -> SetGuidance
    ("Use \"/vis/set/textLayout\" to set layout:");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("x", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("y", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("font_size", 'd', omitable = true);
  parameter->SetDefaultValue (12);
  parameter->SetGuidance ("pixels");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("x_offset", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  parameter->SetGuidance ("pixels");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("y_offset", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  parameter->SetGuidance ("pixels");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("text", 's', omitable = true);
  parameter->SetGuidance ("The rest of the line is text.");
  parameter->SetDefaultValue ("Hello G4");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSceneAddText2D::~G4VisCommandSceneAddText2D () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddText2D::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddText2D::SetNewValue (G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn = verbosity >= G4VisManager::warnings;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4Tokenizer next(newValue);
  G4double x = StoD(next());
  G4double y = StoD(next());
  G4double font_size = StoD(next());
  G4double x_offset = StoD(next());
  G4double y_offset = StoD(next());
  G4String text = next("\n");

  G4Text g4text(text, G4Point3D(x,y,0.));
  G4VisAttributes visAtts(fCurrentTextColour);
  g4text.SetVisAttributes(visAtts);
  g4text.SetLayout(fCurrentTextLayout);
  g4text.SetScreenSize(font_size);
  g4text.SetOffset(x_offset,y_offset);
  G4Text2D* g4text2D = new G4Text2D(g4text);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddText2D::G4Text2D>(g4text2D);
  model->SetType("Text2D");
  model->SetGlobalTag("Text2D");
  model->SetGlobalDescription("Text2D: " + newValue);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "2D text \"" << text
	     << "\" has been added to scene \"" << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

G4VisCommandSceneAddText2D::G4Text2D::G4Text2D(const G4Text& text):
  fText(text)
{}

void G4VisCommandSceneAddText2D::G4Text2D::operator()
  (G4VGraphicsScene& sceneHandler, const G4ModelingParameters*) {
  sceneHandler.BeginPrimitives2D();
  sceneHandler.AddPrimitive(fText);
  sceneHandler.EndPrimitives2D();
}


////////////// /vis/scene/add/trajectories ///////////////////////////////////

G4VisCommandSceneAddTrajectories::G4VisCommandSceneAddTrajectories () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString
    ("/vis/scene/add/trajectories", this);
  fpCommand -> SetGuidance
    ("Adds trajectories to current scene.");
  fpCommand -> SetGuidance
    ("Causes trajectories, if any, to be drawn at the end of processing an"
     "\nevent.  Switches on trajectory storing and sets the"
     "\ndefault trajectory type.");
  fpCommand -> SetGuidance
    ("The command line parameter list determines the default trajectory type."
     "\nIf it contains the string \"smooth\", auxiliary inter-step points will"
     "\nbe inserted to improve the smoothness of the drawing of a curved"
     "\ntrajectory."
     "\nIf it contains the string \"rich\", significant extra information will"
     "\nbe stored in the trajectory (G4RichTrajectory) amenable to modeling"
     "\nand filtering with \"/vis/modeling/trajectories/create/drawByAttribute\""
     "\nand \"/vis/filtering/trajectories/create/attributeFilter\" commands."
     "\nIt may contain both strings in any order.");
  fpCommand -> SetGuidance
    ("\nTo switch off trajectory storing: \"/tracking/storeTrajectory 0\"."
     "\nSee also \"/vis/scene/endOfEventAction\".");
  fpCommand -> SetGuidance
    ("Note:  This only sets the default.  Independently of the result of this"
     "\ncommand, a user may instantiate a trajectory that overrides this default"
     "\nin PreUserTrackingAction.");
  fpCommand -> SetParameterName ("default-trajectory-type", omitable = true);
  fpCommand -> SetDefaultValue ("");
}

G4VisCommandSceneAddTrajectories::~G4VisCommandSceneAddTrajectories () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddTrajectories::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddTrajectories::SetNewValue (G4UIcommand*,
						    G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn = verbosity >= G4VisManager::warnings;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }
  const G4String& currentSceneName = pScene -> GetName ();

  G4bool smooth = false;
  G4bool rich = false;
  if (newValue.find("smooth") != std::string::npos) smooth = true;
  if (newValue.find("rich") != std::string::npos) rich = true;
  if (newValue.size() && !(rich || smooth)) {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: Unrecognised parameter \"" << newValue << "\""
      "\n  No action taken."
      << G4endl;
    }
    return;
  }

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4String defaultTrajectoryType;
  if (smooth && rich) {
    UImanager->ApplyCommand("/tracking/storeTrajectory 4");
    defaultTrajectoryType = "G4RichTrajectory configured for smooth steps";
  } else if (smooth) {
    UImanager->ApplyCommand("/tracking/storeTrajectory 2");
    defaultTrajectoryType = "G4SmoothTrajectory";
  } else if (rich) {
    UImanager->ApplyCommand("/tracking/storeTrajectory 3");
    defaultTrajectoryType = "G4RichTrajectory";
  } else {
    UImanager->ApplyCommand("/tracking/storeTrajectory 1");
    defaultTrajectoryType = "G4Trajectory";
  }

  if (verbosity >= G4VisManager::errors) {
    G4warn <<
      "Attributes available for modeling and filtering with"
      "\n  \"/vis/modeling/trajectories/create/drawByAttribute\" and"
      "\n  \"/vis/filtering/trajectories/create/attributeFilter\" commands:"
	   << G4endl;
    G4warn << *G4TrajectoriesModel().GetAttDefs();
    if (rich) {
      G4warn << *G4RichTrajectory().GetAttDefs()
	     << *G4RichTrajectoryPoint().GetAttDefs();
    } else if (smooth) {
      G4warn << *G4SmoothTrajectory().GetAttDefs()
	     << *G4SmoothTrajectoryPoint().GetAttDefs();
    } else {
      G4warn << *G4Trajectory().GetAttDefs()
	     << *G4TrajectoryPoint().GetAttDefs();
    }
  }

  const auto& eoeList = pScene->GetEndOfEventModelList();
  auto eoeModel = eoeList.begin();
  for (; eoeModel != eoeList.end(); ++eoeModel) {
    const auto* actualModel = eoeModel->fpModel;
    if (dynamic_cast<const G4TrajectoriesModel*>(actualModel)) break;
  }
  if (eoeModel == eoeList.end()) {
    // No trajectories model exists in the scene so create a new one...
    G4VModel* model = new G4TrajectoriesModel();
    pScene -> AddEndOfEventModel (model, warn);
  }  // ...else it already exists and there is no need to add a new one
  // because G4TrajectoriesModel simply describes trajectories in the
  // trajectories store whatever the type.

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Default trajectory type " << defaultTrajectoryType
	   << "\n  will be used to store trajectories for scene \""
	   << currentSceneName << "\"."
	   << G4endl;
  }

  if (verbosity >= G4VisManager::warnings) {
    G4warn <<
      "WARNING: Trajectory storing has been requested.  This action may be"
      "\n  reversed with \"/tracking/storeTrajectory 0\"."
	   << G4endl;
  }

  CheckSceneAndNotifyHandlers (pScene);
}

////////////// /vis/scene/add/userAction ///////////////////////////////////

G4VisCommandSceneAddUserAction::G4VisCommandSceneAddUserAction () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/scene/add/userAction",this);
  fpCommand -> SetGuidance
    ("Add named Vis User Action to current scene.");
  fpCommand -> SetGuidance
    ("Attempts to match search string to name of action - use unique sub-string.");
  fpCommand -> SetGuidance
    ("(Use /vis/list to see names of registered actions.)");
  fpCommand -> SetGuidance
    ("If name == \"all\" (default), all actions are added.");
  fpCommand -> SetParameterName("action-name", omitable = true);
  fpCommand -> SetDefaultValue("all");
}

G4VisCommandSceneAddUserAction::~G4VisCommandSceneAddUserAction () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddUserAction::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddUserAction::SetNewValue
(G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4bool any = false;

  const std::vector<G4VisManager::UserVisAction>& runDurationUserVisActions =
    fpVisManager->GetRunDurationUserVisActions();
  for (size_t i = 0; i < runDurationUserVisActions.size(); i++) {
    const G4String& name = runDurationUserVisActions[i].fName;
    G4VUserVisAction* visAction = runDurationUserVisActions[i].fpUserVisAction;
    if (newValue == "all" || name.find(newValue) != std::string::npos) {
      any = true;
      AddVisAction(name,visAction,pScene,runDuration,verbosity);
    }
  }

  const std::vector<G4VisManager::UserVisAction>& endOfEventUserVisActions =
    fpVisManager->GetEndOfEventUserVisActions();
  for (size_t i = 0; i < endOfEventUserVisActions.size(); i++) {
    const G4String& name = endOfEventUserVisActions[i].fName;
    G4VUserVisAction* visAction = endOfEventUserVisActions[i].fpUserVisAction;
    if (newValue == "all" || name.find(newValue) != std::string::npos) {
      any = true;
      AddVisAction(name,visAction,pScene,endOfEvent,verbosity);
    }
  }

  const std::vector<G4VisManager::UserVisAction>& endOfRunUserVisActions =
    fpVisManager->GetEndOfRunUserVisActions();
  for (size_t i = 0; i < endOfRunUserVisActions.size(); i++) {
    const G4String& name = endOfRunUserVisActions[i].fName;
    G4VUserVisAction* visAction = endOfRunUserVisActions[i].fpUserVisAction;
    if (newValue == "all" || name.find(newValue) != std::string::npos) {
      any = true;
      AddVisAction(name,visAction,pScene,endOfRun,verbosity);
    }
  }

  if (!any) {
    if (verbosity >= G4VisManager::warnings) {
      G4warn <<	"WARNING: No User Vis Action registered." << G4endl;
    }
    return;
  }

  CheckSceneAndNotifyHandlers (pScene);
}

void G4VisCommandSceneAddUserAction::AddVisAction
(const G4String& name,
 G4VUserVisAction* visAction,
 G4Scene* pScene,
 G4VisCommandSceneAddUserAction::ActionType type,
 G4VisManager::Verbosity verbosity)
{
  G4bool warn = verbosity >= G4VisManager::warnings;

  const std::map<G4VUserVisAction*,G4VisExtent>& visExtentMap =
    fpVisManager->GetUserVisActionExtents();
  G4VisExtent extent;
  std::map<G4VUserVisAction*,G4VisExtent>::const_iterator i =
    visExtentMap.find(visAction);
  if (i != visExtentMap.end()) extent = i->second;
  if (warn) {
    if (extent.GetExtentRadius() <= 0.) {
      G4warn
      << "WARNING: User Vis Action \"" << name << "\" extent is null."
      << G4endl;
    }
  }

  G4VModel* model = new G4CallbackModel<G4VUserVisAction>(visAction);
  model->SetType("User Vis Action");
  model->SetGlobalTag(name);
  model->SetGlobalDescription(name);
  model->SetExtent(extent);
  G4bool successful = false;;
  switch (type) {
  case runDuration:
    successful = pScene -> AddRunDurationModel (model, warn);
    break;
  case endOfEvent:
    successful = pScene -> AddEndOfEventModel (model, warn);
    break;
  case endOfRun:
    successful = pScene -> AddEndOfRunModel (model, warn);
    break;
  }
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      const G4String& currentSceneName = pScene -> GetName ();
      G4cout << "User Vis Action added to scene \""
      << currentSceneName << "\"";
      if (verbosity >= G4VisManager::parameters) {
        G4cout << "\n  with extent " << extent;
      }
      G4cout << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
}

////////////// /vis/scene/add/volume ///////////////////////////////////////

G4VisCommandSceneAddVolume::G4VisCommandSceneAddVolume () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/volume", this);
  fpCommand -> SetGuidance 
   ("Adds a physical volume to current scene, with optional clipping volume.");
  fpCommand -> SetGuidance 
    ("If physical-volume-name is \"world\" (the default), the top of the"
     "\nmain geometry tree (material world) is added. If \"worlds\", the"
     "\ntops of all worlds - material world and parallel worlds, if any - are"
     "\nadded. Otherwise a search of all worlds is made.");
  fpCommand -> SetGuidance
    ("In the last case the names of all volumes in all worlds are matched"
     "\nagainst physical-volume-name. If this is of the form \"/regexp/\","
     "\nwhere regexp is a regular expression (see C++ regex), the match uses"
     "\nthe usual rules of regular expression matching. Otherwise an exact"
     "\nmatch is required."
     "\nFor example, \"/Shap/\" adds \"Shape1\" and \"Shape2\".");
  fpCommand -> SetGuidance
    ("It may help to see a textual representation of the geometry hierarchy of"
     "\nthe worlds. Try \"/vis/drawTree [worlds]\" or one of the driver/browser"
     "\ncombinations that have the required functionality, e.g., HepRepFile.");
  fpCommand -> SetGuidance
    ("If clip-volume-type is specified, the subsequent parameters are used to"
     "\nto define a clipping volume. For example,"
     "\n\"/vis/scene/add/volume ! ! ! -box km 0 1 0 1 0 1\" will draw the world"
     "\nwith the positive octant cut away. (If the Boolean Processor issues"
     "\nwarnings try replacing 0 by 0.000000001 or something.)");
  fpCommand -> SetGuidance
    ("If clip-volume-type is prepended with '-', the clip-volume is subtracted"
     "\n(cutaway). (This is the default if there is no prepended character.)"
     "\nIf '*' is prepended, the intersection of the physical-volume and the"
     "\nclip-volume is made. (You can make a section through the detector with"
     "\na thin box, for example).");
  fpCommand -> SetGuidance
    ("For \"box\", the parameters are xmin,xmax,ymin,ymax,zmin,zmax."
     "\nOnly \"box\" is programmed at present.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("physical-volume-name", 's', omitable = true);
  parameter -> SetDefaultValue ("world");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("copy-no", 'i', omitable = true);
  parameter -> SetGuidance ("If negative, matches any copy no.");
  parameter -> SetDefaultValue (-1);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("depth-of-descent", 'i', omitable = true);
  parameter -> SetGuidance
    ("Depth of descent of geometry hierarchy. Default = unlimited depth.");
  parameter -> SetDefaultValue (G4PhysicalVolumeModel::UNLIMITED);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("clip-volume-type", 's', omitable = true);
  parameter -> SetParameterCandidates("none box -box *box");
  parameter -> SetDefaultValue ("none");
  parameter -> SetGuidance("[-|*]type.  See general guidance.");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("parameter-unit", 's', omitable = true);
  parameter -> SetDefaultValue ("m");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("parameter-1", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("parameter-2", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("parameter-3", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("parameter-4", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("parameter-5", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("parameter-6", 'd', omitable = true);
  parameter -> SetDefaultValue (0.);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddVolume::~G4VisCommandSceneAddVolume () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddVolume::GetCurrentValue (G4UIcommand*) {
  return "world 0 -1";
}

void G4VisCommandSceneAddVolume::SetNewValue (G4UIcommand*,
					      G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn = verbosity >= G4VisManager::warnings;

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String name, clipVolumeType, parameterUnit;
  G4int copyNo, requestedDepthOfDescent;
  G4double param1, param2, param3, param4, param5, param6;
  std::istringstream is (newValue);
  is >> name >> copyNo >> requestedDepthOfDescent
     >> clipVolumeType >> parameterUnit
     >> param1 >> param2 >> param3 >> param4 >> param5 >> param6;
  G4PhysicalVolumeModel::ClippingMode clippingMode =
    G4PhysicalVolumeModel::subtraction;  // Default subtraction mode.
  if (clipVolumeType[size_t(0)] == '-') {
    clipVolumeType = clipVolumeType.substr(1);  // Remove first character.
  } else if (clipVolumeType[size_t(0)] == '*') {
    clippingMode = G4PhysicalVolumeModel::intersection;
    clipVolumeType = clipVolumeType.substr(1);
  }
  G4double unit = G4UIcommand::ValueOf(parameterUnit);
  param1 *= unit; param2 *= unit; param3 *= unit;
  param4 *= unit; param5 *= unit; param6 *= unit;

  G4VSolid* clippingSolid = nullptr;
  if (clipVolumeType == "box") {
    const G4double dX = (param2 - param1) / 2.;
    const G4double dY = (param4 - param3) / 2.;
    const G4double dZ = (param6 - param5) / 2.;
    const G4double x0 = (param2 + param1) / 2.;
    const G4double y0 = (param4 + param3) / 2.;
    const G4double z0 = (param6 + param5) / 2.;
    clippingSolid = new G4DisplacedSolid
    ("_displaced_clipping_box",
     new G4Box("_clipping_box",dX,dY,dZ),
     G4Translate3D(x0,y0,z0));
  }

  G4TransportationManager* transportationManager =
    G4TransportationManager::GetTransportationManager ();

  size_t nWorlds = transportationManager->GetNoWorlds();
  if (nWorlds > 1) {  // Parallel worlds in operation...
    if (verbosity >= G4VisManager::warnings) {
      static G4bool warned = false;
      if (!warned && name != "worlds") {
	G4warn <<
	  "WARNING: Parallel worlds in operation.  To visualise, specify"
	  "\n  \"worlds\" or the parallel world volume or sub-volume name"
	  "\n   and control visibility with /vis/geometry."
	       << G4endl;
	std::vector<G4VPhysicalVolume*>::iterator iterWorld =
	  transportationManager->GetWorldsIterator();
	for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
	  G4warn << "  World " << i << ": " << (*iterWorld)->GetName()
		 << G4endl;
	  warned = true;
	}
      }
    }
  }

  // Get the world (the initial value of the iterator points to the mass world).
  G4VPhysicalVolume* world = *(transportationManager->GetWorldsIterator());

  if (!world) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
	"ERROR: G4VisCommandSceneAddVolume::SetNewValue:"
	"\n  No world.  Maybe the geometry has not yet been defined."
	"\n  Try \"/run/initialize\""
	     << G4endl;
    }
    return;
  }

  std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVector;

  // When it comes to determining the extent of a physical volume we normally
  // assume the user wishes to ignore "invisible" volumes. For example, most
  // users make the world volume invisible. So we ask the physical volume
  // model to traverse the geometry hierarchy, starting at the named physical
  // volume, until it finds non-invisible ones, whose extents are accumulated
  // to determine the overall extent. (Once a non-invisible volume is found,
  // the search is curtailed - daughters are always contained within the mother
  // so they have no subsequent influence on the extent of the mother - but the
  // search continues at the same level until all highest level non-invisible
  // volumes are found an their extents accumulated.) So the default is
  G4bool useFullExtent = false;
  // However, the above procedure can be time consuming in some situations, such
  // as a nested parameterisation whose ultimate volumes are the first non-
  // visible ones, which are typical of a medical "phantom". So we assume here
  // below that if a user specifies a name other than "world" or "worlds" he/she
  // wished the extent to be determined by the volume, whether it is visible
  // or not. So we set useFullExtent true at that point below.

  if (name == "world") {

    findingsVector.push_back
    (G4PhysicalVolumesSearchScene::Findings(world,world));

  } else if (name == "worlds") {

    if (nWorlds <= 1) {
      if (verbosity >= G4VisManager::warnings) {
	G4warn <<
	  "WARNING: G4VisCommandSceneAddVolume::SetNewValue:"
	  "\n  Parallel worlds requested but none exist."
	  "\n  Just adding material world."
	       << G4endl;
      }
    }
    std::vector<G4VPhysicalVolume*>::iterator iterWorld =
      transportationManager->GetWorldsIterator();
    for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
      findingsVector.push_back
      (G4PhysicalVolumesSearchScene::Findings
       (*iterWorld,*iterWorld));
    }

  } else {  // Search all worlds...

    // Use the model's full extent. This assumes the user wants these
    // volumes in the findings vector (there could be more than one) to
    // determine the scene's extent. Otherwise G4PhysicalVolumeModel would
    // re-calculate each volume's extent based on visibility, etc., which
    // could be time consuming.
    useFullExtent = true;

    std::vector<G4VPhysicalVolume*>::iterator iterWorld =
      transportationManager->GetWorldsIterator();
    for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
      G4ModelingParameters mp;  // Default - no culling.
      G4PhysicalVolumeModel searchModel
      (*iterWorld,
       requestedDepthOfDescent,
       G4Transform3D(),
       &mp,
       useFullExtent);
      G4PhysicalVolumesSearchScene searchScene(&searchModel, name, copyNo);
      searchModel.DescribeYourselfTo (searchScene);  // Initiate search.
      for (const auto& findings: searchScene.GetFindings()) {
        findingsVector.push_back(findings);
      }
    }
  }

  for (const auto& findings: findingsVector) {
    // Set copy number from search findings for replicas and parameterisations.
    findings.fpFoundPV->SetCopyNo(findings.fFoundPVCopyNo);
    G4PhysicalVolumeModel* foundPVModel = new G4PhysicalVolumeModel
    (findings.fpFoundPV,
     requestedDepthOfDescent,
     findings.fFoundObjectTransformation,
     0, // No modelling parameters (these are set later by the scene handler).
     useFullExtent,
     findings.fFoundBasePVPath);
    if (clippingSolid) {
      foundPVModel->SetClippingSolid(clippingSolid);
      foundPVModel->SetClippingMode(clippingMode);
    }
    if (!foundPVModel->Validate(warn)) return;
    // ...so add it to the scene.
    G4bool successful = pScene->AddRunDurationModel(foundPVModel,warn);
    if (successful) {
      if (verbosity >= G4VisManager::confirmations) {
        G4cout << "\"" << findings.fpFoundPV->GetName()
        << "\", copy no. " << findings.fFoundPVCopyNo
        << ",\n  found in searched volume \""
        << findings.fpSearchPV->GetName()
        << "\" at depth " << findings.fFoundDepth
        << ",\n  base path: \"" << findings.fFoundBasePVPath
        << "\",\n  with a requested depth of further descent of ";
        if (requestedDepthOfDescent < 0) {
          G4cout << "<0 (unlimited)";
        }
        else {
          G4cout << requestedDepthOfDescent;
        }
        G4cout << ",\n  has been added to scene \"" << pScene->GetName() << "\"."
        << G4endl;
      }
    } else {
      G4VisCommandsSceneAddUnsuccessful(verbosity);
    }
  }

  if (findingsVector.empty()) {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: Volume \"" << name << "\"";
      if (copyNo >= 0) {
        G4warn << ", copy no. " << copyNo << ",";
      }
      G4warn << " not found." << G4endl;
    }
    G4VisCommandsSceneAddUnsuccessful(verbosity);
    return;
  }

  CheckSceneAndNotifyHandlers(pScene);
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/scene/add/plotter ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandSceneAddPlotter::G4VisCommandSceneAddPlotter () {
  fpCommand = new G4UIcommand("/vis/scene/add/plotter", this);
  fpCommand -> SetGuidance ("Add a plotter to current scene.");
  
  G4UIparameter* parameter;
  parameter =  new G4UIparameter ("plotter", 's',false);
  fpCommand->SetParameter(parameter);
}

G4VisCommandSceneAddPlotter::~G4VisCommandSceneAddPlotter () {delete fpCommand;}

G4String G4VisCommandSceneAddPlotter::GetCurrentValue (G4UIcommand*) {return "";}

void G4VisCommandSceneAddPlotter::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(newValue);
  G4VModel* model = new G4PlotterModel(_plotter,newValue);

  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddEndOfRunModel(model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout
      << "Plotter \"" << model->GetCurrentDescription()
      << "\" has been added to scene \"" << currentSceneName << "\"."
      << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);

  CheckSceneAndNotifyHandlers (pScene);
}

