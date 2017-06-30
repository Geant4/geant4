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
// $Id: G4VisCommandsSceneAdd.cc 104163 2017-05-15 06:52:42Z gcosmo $
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
#include "G4MagneticFieldModel.hh"
#include "G4PSHitsModel.hh"
#include "G4TrajectoriesModel.hh"
#include "G4ScaleModel.hh"
#include "G4TextModel.hh"
#include "G4ArrowModel.hh"
#include "G4AxesModel.hh"
#include "G4PhysicalVolumeSearchScene.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ApplicationState.hh"
#include "G4VUserVisAction.hh"
#include "G4CallbackModel.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Polyhedron.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4Tokenizer.hh"
#include "G4RunManager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif
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

#include <sstream>

// Local function with some frequently used error printing...
static void G4VisCommandsSceneAddUnsuccessful
(G4VisManager::Verbosity verbosity) {
  if (verbosity >= G4VisManager::warnings) {
    G4cout <<
      "WARNING: For some reason, possibly mentioned above, it has not been"
      "\n  possible to add to the scene."
	   << G4endl;
  }
}

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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
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
  (G4VGraphicsScene& sceneHandler, const G4Transform3D&)
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  } else {
    if (pScene->GetExtent().GetExtentRadius() <= 0.) {
      if (verbosity >= G4VisManager::errors) {
        G4cerr
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
    0.005 * fCurrentLineWidth * sceneExtent.GetExtentRadius();
  // ...but limit it to length/50.
  if (arrowWidth > length/50.) arrowWidth = length/50.;

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
  UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  if (layoutString(0) == 'l') layout = G4Text::left;
  else if (layoutString(0) == 'c') layout = G4Text::centre;
  else if (layoutString(0) == 'r') layout = G4Text::right;

  Date* date = new Date(fpVisManager, size, x, y, layout, dateString);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddDate::Date>(date);
  model->SetType("Date");
  model->SetGlobalTag("Date");
  model->SetGlobalDescription("Date");
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
  UpdateVisManagerScene (currentSceneName);
}

void G4VisCommandSceneAddDate::Date::operator()
  (G4VGraphicsScene& sceneHandler, const G4Transform3D&)
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int size;
  G4double x, y;
  G4String layoutString;
  std::istringstream is(newValue);
  is >> size >> x >> y >> layoutString;

  G4Text::Layout layout = G4Text::right;
  if (layoutString(0) == 'l') layout = G4Text::left;
  else if (layoutString(0) == 'c') layout = G4Text::centre;
  else if (layoutString(0) == 'r') layout = G4Text::right;

  EventID* eventID = new EventID(fpVisManager, size, x, y, layout);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddEventID::EventID>(eventID);
  model->SetType("EventID");
  model->SetGlobalTag("EventID");
  model->SetGlobalDescription("EventID");
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddEndOfEventModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "EventID has been added to scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
  UpdateVisManagerScene (currentSceneName);
}

void G4VisCommandSceneAddEventID::EventID::operator()
  (G4VGraphicsScene& sceneHandler, const G4Transform3D&)
{
  const G4Run* currentRun = 0;
  G4RunManager* runManager = G4RunManager::GetRunManager();
#ifdef G4MULTITHREADED
  if (G4Threading::IsMultithreadedApplication()) {
    runManager = G4MTRunManager::GetMasterRunManager();
  }
#endif
  if (runManager) currentRun = runManager->GetCurrentRun();

  G4VModel* model = fpVisManager->GetCurrentSceneHandler()->GetModel();
  const G4ModelingParameters* mp = 0;
  const G4Event* currentEvent = 0;
  if (model) {
   mp = model->GetModelingParameters();
   currentEvent = mp->GetEvent();
  } else {
    G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<	"ERROR: No model defined for this SceneHandler : "
	     << fpVisManager->GetCurrentSceneHandler()->GetName()
	     << G4endl;
    }
  }
  if (currentRun && currentEvent) {
    G4int runID = currentRun->GetRunID();
    G4int eventID = currentEvent->GetEventID();
    std::ostringstream oss;
    if (fpVisManager->GetCurrentScene()->GetRefreshAtEndOfEvent()) {
      oss << "Run " << runID << " Event " << eventID;
    } else {
      G4int nEvents = 0;
#ifdef G4MULTITHREADED
      if (G4Threading::G4GetThreadId() != G4Threading::MASTER_ID) {
        // Vis sub-thread - run in progress
#else
      G4StateManager* stateManager = G4StateManager::GetStateManager();
      G4ApplicationState state = stateManager->GetCurrentState();
      if (state == G4State_EventProc) {
        // Run in progress
#endif
        nEvents = currentRun->GetNumberOfEventToBeProcessed();
      } else {
        // Rebuilding from kept events
        const std::vector<const G4Event*>* events =
        currentRun->GetEventVector();
        if (events) nEvents = events->size();
      }
#ifndef G4MULTITHREADED
      // In sequential mode we can recognise the last event and avoid
      // drawing the event ID.  But for MT mode there is no way of
      // knowing so we have to accept that the event ID will be drawn
      // at each event, the same characters over and over - but hey!
      if (eventID < nEvents - 1) return;  // Not last event.
#endif
      oss << "Run " << runID << " (" << nEvents << " event";
      if (nEvents != 1) oss << 's';
      oss << ')';
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
}

G4VisCommandSceneAddExtent::Extent::Extent
(G4double xmin, G4double xmax,
 G4double ymin, G4double ymax,
 G4double zmin, G4double zmax):
fExtent(xmin,xmax,ymin,ymax,zmin,zmax)
{}

void G4VisCommandSceneAddExtent::Extent::operator()
(G4VGraphicsScene&, const G4Transform3D&)
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
}

void G4VisCommandSceneAddFrame::Frame::operator()
  (G4VGraphicsScene& sceneHandler, const G4Transform3D&)
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
        G4cerr << "ERROR: No current scene.  Please create one." << G4endl;
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
    UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
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
  (G4VGraphicsScene& sceneHandler, const G4Transform3D&)
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
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
  (G4VGraphicsScene& sceneHandler, const G4Transform3D&)
{
  sceneHandler.BeginPrimitives2D();
  sceneHandler.AddPrimitive(fPolyline);
  sceneHandler.EndPrimitives2D();
}

////////////// /vis/scene/add/logicalVolume //////////////////////////////////

G4VisCommandSceneAddLogicalVolume::G4VisCommandSceneAddLogicalVolume () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/logicalVolume", this);
  fpCommand -> SetGuidance ("Adds a logical volume to the current scene,");
  fpCommand -> SetGuidance
  ("Shows boolean components (if any), voxels (if any), readout geometry"
   "\n  (if any) and local axes, under control of the appropriate flag."
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String name;
  G4int requestedDepthOfDescent;
  G4String booleansString, voxelsString, readoutString, axesString;
  std::istringstream is (newValue);
  is >> name >> requestedDepthOfDescent
     >>  booleansString >> voxelsString >> readoutString >> axesString;
  G4bool booleans = G4UIcommand::ConvertToBool(booleansString);
  G4bool voxels = G4UIcommand::ConvertToBool(voxelsString);
  G4bool readout = G4UIcommand::ConvertToBool(readoutString);
  G4bool axes = G4UIcommand::ConvertToBool(axesString);

  G4LogicalVolumeStore *pLVStore = G4LogicalVolumeStore::GetInstance();
  int nLV = pLVStore -> size ();
  int iLV;
  G4LogicalVolume* pLV = 0;
  for (iLV = 0; iLV < nLV; iLV++ ) {
    pLV = (*pLVStore) [iLV];
    if (pLV -> GetName () == name) break;
  }
  if (iLV == nLV) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Logical volume " << name
	     << " not found in logical volume store." << G4endl;
    }
    return;
  }

  const std::vector<G4Scene::Model>& rdModelList =
    pScene -> GetRunDurationModelList();
  std::vector<G4Scene::Model>::const_iterator i;
  for (i = rdModelList.begin(); i != rdModelList.end(); ++i) {
    if (i->fpModel->GetGlobalDescription().find("Volume") != std::string::npos) break;
  }
  if (i != rdModelList.end()) {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "There is already a volume, \""
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
    (pLV, requestedDepthOfDescent, booleans, voxels, readout);
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
//      G4cout << G4AttCheck(attValues, attDefs);
//      delete attValues;
//    }

    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Logical volume \"" << pLV -> GetName ()
	     << " with requested depth of descent "
	     << requestedDepthOfDescent
	     << ",\n with";
      if (!booleans) G4cout << "out";
      G4cout << " boolean components, with";
      if (!voxels) G4cout << "out";
      G4cout << " voxels and with";
      if (!readout) G4cout << "out";
      G4cout << " readout geometry,"
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

  UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  } else {
    if (pScene->GetExtent().GetExtentRadius() <= 0.) {
      if (verbosity >= G4VisManager::errors) {
        G4cerr
  << "ERROR: Scene has no extent. Add volumes or use \"/vis/scene/add/extent\"."
        << G4endl;
      }
      return;
    }
  }

  G4VViewer* pViewer = fpVisManager->GetCurrentViewer();
  if (!pViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
	"ERROR: G4VisCommandSceneAddLogo::SetNewValue: no viewer."
	"\n  Auto direction needs a viewer."
	     << G4endl;
    }
    return;
  }

  G4double userHeight, red, green, blue, xmid, ymid, zmid;
  G4String userHeightUnit, direction, auto_manual, positionUnit;
  std::istringstream is (newValue);
  is >> userHeight >> userHeightUnit >> direction
     >> red >> green >> blue
     >> auto_manual
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
  else if (direction(0) == 'x') logoDirection = X;
  else if (direction(0) == 'y') logoDirection = Y;
  else if (direction(0) == 'z') logoDirection = Z;
  else if (direction(0) == '-') {
    if (direction(1) == 'x') logoDirection = minusX;
    else if (direction(1) == 'y') logoDirection = minusY;
    else if (direction(1) == 'z') logoDirection = minusZ;
  } else {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<	"ERROR: Unrecogniseed direction: \""
	     << direction << "\"." << G4endl;
      return;
    }
  }

  G4bool autoPlacing = false; if (auto_manual == "auto") autoPlacing = true;
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
      G4cout <<
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
      G4cout <<
	"WARNING: Not enough room in existing scene.  Maybe logo is too large."
	     << G4endl;
    }
  }
  if (worried) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
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

  G4Logo* logo = new G4Logo(height,visAtts);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddLogo::G4Logo>(logo);
  model->SetType("G4Logo");
  model->SetGlobalTag("G4Logo");
  model->SetGlobalDescription("G4Logo: " + newValue);
  model->SetTransformation(transform);
  // Note: it is the responsibility of the model to act upon this, but
  // the extent is in local coordinates...
  G4double& h = height;
  G4double h2 = h/2.;
  G4VisExtent extent(-h,h,-h2,h2,-h2,h2);
  model->SetExtent(extent);
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
	       << transform.getTranslation();
      }
      G4cout << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
  UpdateVisManagerScene (currentSceneName);
}

G4VisCommandSceneAddLogo::G4Logo::G4Logo
(G4double height, const G4VisAttributes& visAtts):
  fVisAtts(visAtts)
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
  fpG->SetVisAttributes(&fVisAtts);
  fpG->Transform(G4Translate3D(-0.55*h,0.,0.));

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
  fp4->SetVisAttributes(&fVisAtts);
  fp4->Transform(G4Translate3D(0.55*h,0.,0.));
}

G4VisCommandSceneAddLogo::G4Logo::~G4Logo() {
  delete fpG;
  delete fp4;
}

void G4VisCommandSceneAddLogo::G4Logo::operator()
  (G4VGraphicsScene& sceneHandler, const G4Transform3D& transform) {
  sceneHandler.BeginPrimitives(transform);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int size;
  G4double x, y;
  G4String layoutString;
  std::istringstream is(newValue);
  is >> size >> x >> y >> layoutString;
  G4Text::Layout layout = G4Text::right;
  if (layoutString(0) == 'l') layout = G4Text::left;
  else if (layoutString(0) == 'c') layout = G4Text::centre;
  else if (layoutString(0) == 'r') layout = G4Text::right;

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
  UpdateVisManagerScene (currentSceneName);
}

void G4VisCommandSceneAddLogo2D::Logo2D::operator()
  (G4VGraphicsScene& sceneHandler, const G4Transform3D&)
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
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/magneticField", this);
  fpCommand -> SetGuidance
  ("Adds magnetic field representation to current scene.");
  fpCommand -> SetGuidance
  ("The first parameter is no. of data points per half scene.  So, possibly, at"
   "\nmaximum, the number of data points sampled is (2*n+1)^3, which can grow"
   "\nlarge--be warned!"
   "\nYou might find that your scene is cluttered by thousands of arrows for"
   "\nthe default number of data points, so try reducing to 2 or 3, e.g:"
   "\n  /vis/scene/add/magneticField 3"
   "\nor, if only a small part of the scene has a field:"
   "\n  /vis/scene/add/magneticField 50 or more");
  fpCommand -> SetGuidance
  ("In the arrow representation, the length of the arrow is proportional"
   "\nto the magnitude of the field and the colour is mapped onto the range"
   "\nas a fraction of the maximum magnitude: 0->0.5->1 is blue->green->red.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("nDataPointsPerHalfScene", 'i', omitable = true);
  parameter -> SetDefaultValue (10);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("representation", 's', omitable = true);
  parameter -> SetParameterCandidates("fullArrow lightArrow");
  parameter -> SetDefaultValue ("fullArrow");
  fpCommand -> SetParameter (parameter);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int nDataPointsPerHalfScene;
  G4String representation;
  std::istringstream iss(newValue);
  iss >> nDataPointsPerHalfScene >> representation;
  G4MagneticFieldModel::Representation
  modelRepresentation = G4MagneticFieldModel::fullArrow;
  if (representation == "lightArrow") {
    modelRepresentation = G4MagneticFieldModel::lightArrow;
  }
  G4VModel* model =
  new G4MagneticFieldModel(nDataPointsPerHalfScene,modelRepresentation,
                           fCurrentArrow3DLineSegmentsPerCircle);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Magnetic field, if any, will be drawn in scene \""
      << currentSceneName
      << "\"\n  with "
      << nDataPointsPerHalfScene
      << " data points per half scene and with representation \""
      << representation
      << '\"'
      << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
  UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
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
  fpCommand -> SetGuidance (G4Scale::GetGuidanceString());
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  } else {
    if (pScene->GetExtent().GetExtentRadius() <= 0.) {
      if (verbosity >= G4VisManager::errors) {
        G4cerr
  << "ERROR: Scene has no extent. Add volumes or use \"/vis/scene/add/extent\"."
        << G4endl;
      }
      return;
    }
  }

  G4double userLength, red, green, blue, xmid, ymid, zmid;
  G4String userLengthUnit, direction, auto_manual, positionUnit;
  std::istringstream is (newValue);
  is >> userLength >> userLengthUnit >> direction
     >> red >> green >> blue
     >> auto_manual
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

  G4Scale::Direction scaleDirection (G4Scale::x);
  if (direction(0) == 'y') scaleDirection = G4Scale::y;
  if (direction(0) == 'z') scaleDirection = G4Scale::z;

  G4VViewer* pViewer = fpVisManager->GetCurrentViewer();
  if (!pViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
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
      if (std::abs(up.y()) > std::abs(up.z())) scaleDirection = G4Scale::z;
	  else scaleDirection = G4Scale::y;
    }
    else if (std::abs(vp.y()) > std::abs(vp.x()) &&
	     std::abs(vp.y()) > std::abs(vp.z())) {  // y viewpoint
      if (std::abs(up.x()) > std::abs(up.z())) scaleDirection = G4Scale::z;
	  else scaleDirection = G4Scale::x;
    }
    else if (std::abs(vp.z()) > std::abs(vp.x()) &&
	     std::abs(vp.z()) > std::abs(vp.y())) {  // z viewpoint
      if (std::abs(up.y()) > std::abs(up.x())) scaleDirection = G4Scale::x;
	  else scaleDirection = G4Scale::y;
    }
  }

  G4bool autoPlacing = false; if (auto_manual == "auto") autoPlacing = true;
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
      G4cout <<
	"WARNING: Existing scene does not yet have any extent."
	"\n  Maybe you have not yet added any geometrical object."
	     << G4endl;
    }
  }
  // Test existing scene for room...
  G4bool room  = true;
  switch (scaleDirection) {
  case G4Scale::x:
    if (freeLengthFraction * (xmax - xmin) < length) room = false;
    break;
  case G4Scale::y:
    if (freeLengthFraction * (ymax - ymin) < length) room = false;
    break;
  case G4Scale::z:
    if (freeLengthFraction * (zmax - zmin) < length) room = false;
    break;
  }
  if (!room) {
    worried = true;
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"WARNING: Not enough room in existing scene.  Maybe scale is too long."
	     << G4endl;
    }
  }
  if (worried) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"WARNING: The scale you have asked for is bigger than the existing"
	"\n  scene.  Maybe you have added it too soon.  It is recommended that"
	"\n  you add the scale last so that it can be correctly auto-positioned"
	"\n  so as not to be obscured by any existing object and so that the"
	"\n  view parameters can be correctly recalculated."
	     << G4endl;
    }
  }

  // Let's go ahead a construct a scale and a scale model.  Since the
  // placing is done here, this G4Scale is *not* auto-placed...
  G4Scale scale(length, annotation, scaleDirection,
		false, xmid, ymid, zmid,
                fCurrentTextSize);
  G4VisAttributes visAttr(G4Colour(red, green, blue));
  scale.SetVisAttributes(visAttr);
  G4VModel* model = new G4ScaleModel(scale);
  G4String globalDescription = model->GetGlobalDescription();
  globalDescription += " (" + newValue + ")";
  model->SetGlobalDescription(globalDescription);

  // Now figure out the extent...
  //
  // From the G4Scale.hh:
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
  // End of clip from G4Scale.hh:
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
    case G4Scale::x:
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
    case G4Scale::y:
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
    case G4Scale::z:
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

  /* Old code - kept for future reference.
  G4double sxmid(xmid), symid(ymid), szmid(zmid);
  if (autoPlacing) {
    sxmid = xmin + onePlusComfort * (xmax - xmin);
    symid = ymin - comfort * (ymax - ymin);
    szmid = zmin + onePlusComfort * (zmax - zmin);
    switch (scaleDirection) {
    case G4Scale::x:
      sxmid -= halfLength;
      break;
    case G4Scale::y:
      symid += halfLength;
      break;
    case G4Scale::z:
      szmid -= halfLength;
      break;
    }
  }
  */

  /* sxmin, etc., not actually used.  Comment out to prevent compiler
     warnings but keep in case need in future.  Extract transform and
     scaleExtent into reduced code below.
  G4double sxmin(sxmid), sxmax(sxmid);
  G4double symin(symid), symax(symid);
  G4double szmin(szmid), szmax(szmid);
  G4Transform3D transform;
  G4VisExtent scaleExtent;
  switch (scaleDirection) {
  case G4Scale::x:
    sxmin = sxmid - halfLength;
    sxmax = sxmid + halfLength;
    scaleExtent = G4VisExtent(-halfLength,halfLength,0,0,0,0);
    break;
  case G4Scale::y:
    symin = symid - halfLength;
    symax = symid + halfLength;
    transform = G4RotateZ3D(halfpi);
    scaleExtent = G4VisExtent(0,0,-halfLength,halfLength,0,0);
    break;
  case G4Scale::z:
    szmin = szmid - halfLength;
    szmax = szmid + halfLength;
    transform = G4RotateY3D(halfpi);
    scaleExtent = G4VisExtent(0,0,0,0,-halfLength,halfLength);
    break;
  }
  */
  G4Transform3D transform;
  G4VisExtent scaleExtent;
  switch (scaleDirection) {
  case G4Scale::x:
    scaleExtent = G4VisExtent(-halfLength,halfLength,0,0,0,0);
    break;
  case G4Scale::y:
    transform = G4RotateZ3D(halfpi);
    scaleExtent = G4VisExtent(0,0,-halfLength,halfLength,0,0);
    break;
  case G4Scale::z:
    transform = G4RotateY3D(halfpi);
    scaleExtent = G4VisExtent(0,0,0,0,-halfLength,halfLength);
    break;
  }
  transform = G4Translate3D(sxmid,symid,szmid) * transform;
  /////////  G4VisExtent scaleExtent(sxmin, sxmax, symin, symax, szmin, szmax);


  model->SetTransformation(transform);
  // Note: it is the responsibility of the model to act upon this, but
  // the extent is in local coordinates...
  model->SetExtent(scaleExtent);
  // This extent gets "added" to existing scene extent in
  // AddRunDurationModel below.

  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Scale of " << annotation
	     << " added to scene \"" << currentSceneName << "\".";
      if (verbosity >= G4VisManager::parameters) {
	G4cout << "\n  with extent " << scaleExtent
	       << "\n  at " << transform.getRotation()
	       << transform.getTranslation();
      }
      G4cout << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
  UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
  UpdateVisManagerScene (currentSceneName);
}

G4VisCommandSceneAddText2D::G4Text2D::G4Text2D(const G4Text& text):
  fText(text)
{}

void G4VisCommandSceneAddText2D::G4Text2D::operator()
  (G4VGraphicsScene& sceneHandler, const G4Transform3D& transform) {
  sceneHandler.BeginPrimitives2D(transform);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4bool smooth = false;
  G4bool rich = false;
  if (newValue.find("smooth") != std::string::npos) smooth = true;
  if (newValue.find("rich") != std::string::npos) rich = true;
  if (newValue.size() && !(rich || smooth)) {
    if (verbosity >= G4VisManager::errors) {
      G4cerr << "ERROR: Unrecognised parameter \"" << newValue << "\""
      "\n  No action taken."
      << G4endl;
    }
    return;
  }

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
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
  UImanager->SetVerboseLevel(keepVerbose);

  if (verbosity >= G4VisManager::errors) {
    G4cout <<
      "Attributes available for modeling and filtering with"
      "\n  \"/vis/modeling/trajectories/create/drawByAttribute\" and"
      "\n  \"/vis/filtering/trajectories/create/attributeFilter\" commands:"
	   << G4endl;
    G4cout << *G4TrajectoriesModel().GetAttDefs();
    if (rich) {
      G4cout << *G4RichTrajectory().GetAttDefs()
	     << *G4RichTrajectoryPoint().GetAttDefs();
    } else if (smooth) {
      G4cout << *G4SmoothTrajectory().GetAttDefs()
	     << *G4SmoothTrajectoryPoint().GetAttDefs();
    } else {
      G4cout << *G4Trajectory().GetAttDefs()
	     << *G4TrajectoryPoint().GetAttDefs();
    }
  }

  G4VModel* model = new G4TrajectoriesModel();
  const G4String& currentSceneName = pScene -> GetName ();
  pScene -> AddEndOfEventModel (model, warn);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Default trajectory type " << defaultTrajectoryType
	   << "\n  will be used to store trajectories for scene \""
	   << currentSceneName << "\"."
	   << G4endl;
  }

  if (verbosity >= G4VisManager::warnings) {
    G4cout <<
      "WARNING: Trajectory storing has been requested.  This action may be"
      "\n  reversed with \"/tracking/storeTrajectory 0\"."
	   << G4endl;
  }
  UpdateVisManagerScene (currentSceneName);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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
      G4cout <<	"WARNING: No User Vis Action registered." << G4endl;
    }
    return;
  }

  const G4String& currentSceneName = pScene -> GetName ();
  UpdateVisManagerScene (currentSceneName);
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
      G4cout <<	"WARNING: User Vis Action extent is null." << G4endl;
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
  if (successful && verbosity >= G4VisManager::confirmations) {
    const G4String& currentSceneName = pScene -> GetName ();
    G4cout << "User Vis Action added to scene \""
	   << currentSceneName << "\"";
    if (verbosity >= G4VisManager::parameters) {
      G4cout << "\n  with extent " << extent;
    }
    G4cout << G4endl;
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
     "\nmain geometry tree (material world) is added.  If \"worlds\", the"
     "\ntop of all worlds - material world and parallel worlds, if any - are"
     "\nadded.  Otherwise a search of all worlds is made, taking the first"
     "\nmatching occurence only.  To see a representation of the geometry"
     "\nhierarchy of the worlds, try \"/vis/drawTree [worlds]\" or one of the"
     "\ndriver/browser combinations that have the required functionality, e.g., HepRep.");
  fpCommand -> SetGuidance
    ("If clip-volume-type is specified, the subsequent parameters are used to"
     "\nto define a clipping volume.  For example,"
     "\n\"/vis/scene/add/volume ! ! ! -box km 0 1 0 1 0 1\" will draw the world"
     "\nwith the positive octant cut away.  (If the Boolean Processor issues"
     "\nwarnings try replacing 0 by 0.000000001 or something.)");
  fpCommand -> SetGuidance
    ("If clip-volume-type is prepended with '-', the clip-volume is subtracted"
     "\n(cutaway). (This is the default if there is no prepended character.)"
     "\nIf '*' is prepended, the intersection of the physical-volume and the"
     "\nclip-volume is made. (You can make a section/DCUT with a thin box, for"
     "\nexample).");
  fpCommand -> SetGuidance
    ("For \"box\", the parameters are xmin,xmax,ymin,ymax,zmin,zmax."
     "\nOnly \"box\" is programmed at present.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("physical-volume-name", 's', omitable = true);
  parameter -> SetDefaultValue ("world");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("copy-no", 'i', omitable = true);
  parameter -> SetGuidance
    ("If negative, matches any copy no.  First name match is taken.");
  parameter -> SetDefaultValue (-1);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("depth-of-descent", 'i', omitable = true);
  parameter -> SetGuidance
    ("Depth of descent of geometry hierarchy. Default = unlimited depth.");
  parameter -> SetDefaultValue (G4Scene::UNLIMITED);
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
      G4cerr <<	"ERROR: No current scene.  Please create one." << G4endl;
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

  G4TransportationManager* transportationManager =
    G4TransportationManager::GetTransportationManager ();

  size_t nWorlds = transportationManager->GetNoWorlds();
  if (nWorlds > 1) {  // Parallel worlds in operation...
    if (verbosity >= G4VisManager::warnings) {
      static G4bool warned = false;
      if (!warned && name != "worlds") {
	G4cout <<
	  "WARNING: Parallel worlds in operation.  To visualise, specify"
	  "\n  \"worlds\" or the parallel world volume or sub-volume name"
	  "\n   and control visibility with /vis/geometry."
	       << G4endl;
	std::vector<G4VPhysicalVolume*>::iterator iterWorld =
	  transportationManager->GetWorldsIterator();
	for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
	  G4cout << "  World " << i << ": " << (*iterWorld)->GetName()
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
      G4cerr <<
	"ERROR: G4VisCommandSceneAddVolume::SetNewValue:"
	"\n  No world.  Maybe the geometry has not yet been defined."
	"\n  Try \"/run/initialize\""
	     << G4endl;
    }
    return;
  }

  std::vector<G4PhysicalVolumeModel*> models;
  std::vector<G4VPhysicalVolume*> foundVolumes;
  G4VPhysicalVolume* foundWorld = 0;
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
  PVPath foundFullPVPath;
  std::vector<G4int> foundDepths;
  std::vector<G4Transform3D> transformations;

  if (name == "world") {

    models.push_back
      (new G4PhysicalVolumeModel (world, requestedDepthOfDescent));
    foundVolumes.push_back(world);
    foundDepths.push_back(0);
    transformations.push_back(G4Transform3D());

  } else if (name == "worlds") {

    if (nWorlds == 0) {
      if (verbosity >= G4VisManager::warnings) {
	G4cout <<
	  "WARNING: G4VisCommandSceneAddVolume::SetNewValue:"
	  "\n  Parallel worlds requested but none exist."
	  "\n  Just adding material world."
	       << G4endl;
      }
    }
    std::vector<G4VPhysicalVolume*>::iterator iterWorld =
      transportationManager->GetWorldsIterator();
    for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
      models.push_back
	(new G4PhysicalVolumeModel (*iterWorld, requestedDepthOfDescent));
      foundVolumes.push_back(*iterWorld);
      foundDepths.push_back(0);
      transformations.push_back(G4Transform3D());
    }

  } else {  // Search all worlds...
    
    std::vector<G4VPhysicalVolume*>::iterator iterWorld =
      transportationManager->GetWorldsIterator();
    for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
      G4PhysicalVolumeModel searchModel (*iterWorld);  // Unlimited depth.
      G4ModelingParameters mp;  // Default - no culling.
      searchModel.SetModelingParameters (&mp);
      G4PhysicalVolumeSearchScene searchScene (&searchModel, name, copyNo);
      searchModel.DescribeYourselfTo (searchScene);  // Initiate search.
      G4VPhysicalVolume* foundVolume = searchScene.GetFoundVolume ();
      if (foundVolume) {
	foundWorld = *iterWorld;
	foundVolumes.push_back(foundVolume);
        foundFullPVPath = searchScene.GetFoundFullPVPath();
        foundDepths.push_back(searchScene.GetFoundDepth());
	transformations.push_back(searchScene.GetFoundTransformation());
	break;
      }
    }

    if (foundVolumes.size()) {
      for (size_t i = 0; i < foundVolumes.size(); ++i) {
        G4PhysicalVolumeModel* foundPVModel = new G4PhysicalVolumeModel
        (foundVolumes[i], requestedDepthOfDescent, transformations[i]);
        foundFullPVPath.pop_back();  // "Base" is "Found - 1".
        foundPVModel->SetBaseFullPVPath(foundFullPVPath);
	models.push_back(foundPVModel);
      }
    } else {
      if (verbosity >= G4VisManager::errors) {
	G4cerr << "ERROR: Volume \"" << name << "\"";
	if (copyNo >= 0) {
	  G4cerr << ", copy no. " << copyNo << ",";
	}
	G4cerr << " not found." << G4endl;
      }
      return;
    }
  }

  if (clipVolumeType == "box") {
    const G4double dX = (param2 - param1) / 2.;
    const G4double dY = (param4 - param3) / 2.;
    const G4double dZ = (param6 - param5) / 2.;
    const G4double x0 = (param2 + param1) / 2.;
    const G4double y0 = (param4 + param3) / 2.;
    const G4double z0 = (param6 + param5) / 2.;
    G4VSolid* clippingSolid =
      new G4DisplacedSolid
      ("_displaced_clipping_box",
       new G4Box("_clipping_box",dX,dY,dZ),
       G4Translate3D(x0,y0,z0));
    for (size_t i = 0; i < foundVolumes.size(); ++i) {
      models[i]->SetClippingSolid(clippingSolid);
      models[i]->SetClippingMode(clippingMode);
    }
  }  // If any other shape consider NumberOfRotationSides!!!!!!!!!!!

  const G4String& currentSceneName = pScene -> GetName ();
  G4bool failure = true;
  for (size_t i = 0; i < foundVolumes.size(); ++i) {
    G4bool successful = pScene -> AddRunDurationModel (models[i], warn);
    if (successful) {
      failure = false;
      if (verbosity >= G4VisManager::confirmations) {
	G4cout << "First occurrence of \""
	       << foundVolumes[i] -> GetName ()
	       << "\"";
	if (copyNo >= 0) {
	  G4cout << ", copy no. " << copyNo << ",";
	}
	G4cout << "\n  found ";
	if (foundWorld)
	  G4cout << "in world \"" << foundWorld->GetName() << "\" ";
	G4cout << "at depth " << foundDepths[i]
	       << ",\n  with a requested depth of further descent of ";
	if (requestedDepthOfDescent < 0) {
	  G4cout << "<0 (unlimited)";
	}
	else {
	  G4cout << requestedDepthOfDescent;
	}
	G4cout << ",\n  has been added to scene \"" << currentSceneName << "\"."
	       << G4endl;
      }
    }
  }

  if (failure) {
    G4VisCommandsSceneAddUnsuccessful(verbosity);
    return;
  }

  UpdateVisManagerScene (currentSceneName);
}
