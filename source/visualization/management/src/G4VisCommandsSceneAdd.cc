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
// $Id: G4VisCommandsSceneAdd.cc,v 1.84 2010-11-06 18:34:26 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// /vis/scene commands - John Allison  9th August 1998

#include "G4VisCommandsSceneAdd.hh"

#include "G4TransportationManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4ModelingParameters.hh"
#include "G4HitsModel.hh"
#include "G4DigiModel.hh"
#include "G4PSHitsModel.hh"
#include "G4TrajectoriesModel.hh"
#include "G4ScaleModel.hh"
#include "G4TextModel.hh"
#include "G4AxesModel.hh"
#include "G4PhysicalVolumeSearchScene.hh"
#include "G4VGlobalFastSimulationManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4FlavoredParallelWorldModel.hh"
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
#include "G4Tokenizer.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4IdentityTrajectoryFilter.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4RichTrajectory.hh"
#include "G4RichTrajectoryPoint.hh"
#include "G4AttDef.hh"
#include "G4ios.hh"
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

////////////// /vis/scene/add/axes //////////////////////////////////

G4VisCommandSceneAddAxes::G4VisCommandSceneAddAxes () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/axes", this);
  fpCommand -> SetGuidance ("Add axes.");
  fpCommand -> SetGuidance
    ("Draws axes at (x0, y0, z0) of given length.");
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
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue  ("m");
  fpCommand->SetParameter     (parameter);
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String unitString;
  G4double x0, y0, z0, length;
  std::istringstream is (newValue);
  is >> x0 >> y0 >> z0 >> length >> unitString;

  G4double unit = G4UIcommand::ValueOf(unitString);
  x0 *= unit; y0 *= unit; z0 *= unit; length *= unit;

  G4VModel* model = new G4AxesModel(x0, y0, z0, length);

  model->SetExtent(G4VisExtent(x0 - length, x0 + length,
			       y0 - length, y0 + length,
			       z0 - length, z0 + length));
  // This extent gets "added" to existing scene extent in
  // AddRunDurationModel below.

  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Axes have been added to scene \"" << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
  UpdateVisManagerScene (currentSceneName);
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4DigiModel* model = new G4DigiModel;
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddEndOfEventModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Digis will be drawn in scene \""
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int size;
  G4double x, y;
  std::istringstream is(newValue);
  is >> size >> x >> y;

  EventID* eventID = new EventID(fpVisManager, size, x, y);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddEventID::EventID>(eventID);
  model->SetGlobalDescription("EventID");
  model->SetGlobalTag("EventID");
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddEndOfEventModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "EventID will be drawn in scene \""
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
      G4cout <<	"ERROR: No model defined for this SceneHandler : "
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
      G4StateManager* stateManager = G4StateManager::GetStateManager();
      G4ApplicationState state = stateManager->GetCurrentState();
      if (state == G4State_EventProc) {
	nEvents = currentRun->GetNumberOfEventToBeProcessed();
      } else {
	const std::vector<const G4Event*>* events =
	  currentRun->GetEventVector();
	if (events) nEvents = events->size();
      }
      if (eventID < nEvents - 1) return;  // Not last event.
      else {
	oss << "Run " << runID << " (" << nEvents << " accumulated events)";
      }
    }
    G4Text text(oss.str(), G4Point3D(fX, fY, 0.));
    text.SetScreenSize(fSize);
    G4VisAttributes textAtts(G4Colour(0.,1.,1));
    text.SetVisAttributes(textAtts);
    sceneHandler.BeginPrimitives2D();
    sceneHandler.AddPrimitive(text);
    sceneHandler.EndPrimitives2D();
  }
}

////////////// /vis/scene/add/ghosts ///////////////////////////////////////

G4VisCommandSceneAddGhosts::G4VisCommandSceneAddGhosts () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/add/ghosts", this);
  fpCommand -> SetGuidance
    ("Adds ghost volumes (G4FlavoredParallelWorld) to the current scene.");
  fpCommand -> SetGuidance ("Selects by particle.");
  fpCommand -> SetParameterName ("particle", omitable = true);
  fpCommand -> SetDefaultValue ("all");
}

G4VisCommandSceneAddGhosts::~G4VisCommandSceneAddGhosts () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddGhosts::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddGhosts::SetNewValue(G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }
  const G4String& currentSceneName = pScene -> GetName ();

  // Gets the G4GlobalFastSimulationManager pointer if any.
  G4VGlobalFastSimulationManager* theGlobalFastSimulationManager;
  if(!(theGlobalFastSimulationManager = 
       G4VGlobalFastSimulationManager::GetConcreteInstance ())){
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: no G4GlobalFastSimulationManager" << G4endl;
    }
    return;
  }
  
  // Gets the G4ParticleTable pointer.
  G4ParticleTable* theParticleTable=G4ParticleTable::GetParticleTable();
  
  // If "all" (the default) loops on all known particles
  if(newValue=="all") 
    {
      G4VFlavoredParallelWorld* CurrentFlavoredWorld = 0;
      G4bool successful = false;
      for (G4int iParticle=0; iParticle<theParticleTable->entries(); 
	   iParticle++)
	{
	  CurrentFlavoredWorld = theGlobalFastSimulationManager->
	    GetFlavoredWorldForThis(theParticleTable->GetParticle(iParticle));
	  
	  if(CurrentFlavoredWorld)
	    successful = successful || pScene -> 
	      AddRunDurationModel(new G4FlavoredParallelWorldModel 
				  (CurrentFlavoredWorld), warn);
	}
      if (successful) 
	{
	  if (verbosity >= G4VisManager::confirmations) 
	    G4cout << "Ghosts have been added to scene \""
		   << currentSceneName << "\"."
		   << G4endl;
	  UpdateVisManagerScene (currentSceneName);
	}
      else 
	{
	  G4cout << "ERROR: There are no ghosts."<<G4endl;
	  G4VisCommandsSceneAddUnsuccessful(verbosity);
	}
      return;
    }
  
  // Given a particle name looks just for the concerned Ghosts, if any.
  G4ParticleDefinition* currentParticle = 
    theParticleTable->FindParticle(newValue);
  
  if (currentParticle == NULL) 
    {
      if (verbosity >= G4VisManager::errors) 
	G4cout << "ERROR: \"" << newValue
	       << "\": not found this particle name!" << G4endl;
      return;
    }
  
  G4VFlavoredParallelWorld* worldForThis =
    theGlobalFastSimulationManager->GetFlavoredWorldForThis(currentParticle);
  if(worldForThis) 
    {
      G4bool successful = pScene -> AddRunDurationModel
	(new G4FlavoredParallelWorldModel (worldForThis), warn);
      if (successful) {
	if (verbosity >= G4VisManager::confirmations) 
	  G4cout << "Ghosts have been added to scene \""
		 << currentSceneName << "\"."
		 << G4endl;
	UpdateVisManagerScene (currentSceneName);
      }
    }
  else 
    if (verbosity >= G4VisManager::errors) 
      {
	G4cout << "ERROR: There are no ghosts for \""<<newValue<<"\""<<G4endl;
	G4VisCommandsSceneAddUnsuccessful(verbosity);
      }
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4HitsModel* model = new G4HitsModel;
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddEndOfEventModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Hits will be drawn in scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
  UpdateVisManagerScene (currentSceneName);
}

////////////// /vis/scene/add/logicalVolume //////////////////////////////////

G4VisCommandSceneAddLogicalVolume::G4VisCommandSceneAddLogicalVolume () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/logicalVolume", this);
  fpCommand -> SetGuidance ("Adds a logical volume to the current scene,");
  fpCommand -> SetGuidance
    ("Shows boolean components (if any), voxels (if any) and readout geometry"
     "\n(if any).  Note: voxels are not constructed until start of run -"
     "\n \"/run/beamOn\".");
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String name;
  G4int requestedDepthOfDescent;
  G4String booleansString, voxelsString, readoutString;
  std::istringstream is (newValue);
  is >> name >> requestedDepthOfDescent
     >>  booleansString >> voxelsString >> readoutString;
  G4bool booleans = G4UIcommand::ConvertToBool(booleansString);
  G4bool voxels = G4UIcommand::ConvertToBool(voxelsString);
  G4bool readout = G4UIcommand::ConvertToBool(readoutString);

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
      G4cout << "ERROR: Logical volume " << name
	     << " not found in logical volume store." << G4endl;
    }
    return;
  }

  const std::vector<G4VModel*>& rdModelList = pScene -> GetRunDurationModelList();
  std::vector<G4VModel*>::const_iterator i;
  for (i = rdModelList.begin(); i != rdModelList.end(); ++i) {
    if ((*i) -> GetGlobalDescription().find("Volume") != std::string::npos) break;
  }
  if (i != rdModelList.end()) {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "There is already a volume, \""
             << (*i) -> GetGlobalDescription()
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
  }

  G4VModel* model = new G4LogicalVolumeModel
    (pLV, requestedDepthOfDescent, booleans, voxels, readout);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
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
	     << "\n  has been added to scene \"" << currentSceneName << "\"."
	     << G4endl;
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
  fpCommand -> SetGuidance 
    ("Adds a G4 logo to the current scene.");
  fpCommand -> SetGuidance 
    ("The placement, if automatic, is similar to that of scale -"
     "\n\"help /vis/scene/add/scale\" for more information.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("height", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("direction", 's', omitable = true);
  parameter->SetGuidance ("'x', 'y' or 'z' - otherwise defaults to 'x'.");
  parameter->SetDefaultValue ("x");
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
  parameter =  new G4UIparameter ("auto|manual", 's', omitable = true);
  parameter->SetGuidance
    ("Automatic placement or manual placement at (xmid,ymid,zmid).");
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
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

  G4double height = userHeight * G4UIcommand::ValueOf(userHeightUnit);
  G4double unit = G4UIcommand::ValueOf(positionUnit);
  xmid *= unit; ymid *= unit; zmid *= unit;

  G4Scale::Direction logoDirection (G4Scale::x);
  if (direction(0) == 'y') logoDirection = G4Scale::y;
  if (direction(0) == 'z') logoDirection = G4Scale::z;

  G4bool autoPlacing = false; if (auto_manual == "auto") autoPlacing = true;
  // Parameters read and interpreted.

  // Useful constants, etc...
  const G4double halfHeight(height / 2.);
  const G4double comfort(0.01);  // 0.15 seems too big.  0.05 might be better.
  const G4double onePlusComfort(1. + comfort);
  const G4double freeHeightFraction (1. + 2. * comfort);

  const G4VisExtent& sceneExtent = pScene->GetExtent();  // Existing extent.
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
  G4bool room = true;
  switch (logoDirection) {
  case G4Scale::x:
    if (freeHeightFraction * (xmax - xmin) < height) room = false; break;
  case G4Scale::y:
    if (freeHeightFraction * (ymax - ymin) < height) room = false; break;
  case G4Scale::z:
    if (freeHeightFraction * (zmax - zmin) < height) room = false; break;
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

  G4double sxmid(xmid), symid(ymid), szmid(zmid);
  if (autoPlacing) {
    sxmid = xmin + onePlusComfort * (xmax - xmin);
    symid = ymin - comfort * (ymax - ymin);
    szmid = zmin + onePlusComfort * (zmax - zmin);
    switch (logoDirection) {
    case G4Scale::x:
      sxmid -= halfHeight;
      break;
    case G4Scale::y:
      symid += halfHeight;
      break;
    case G4Scale::z:
      szmid -= halfHeight;
      break;
    }
  }
  G4double sxmin(sxmid), sxmax(sxmid);
  G4double symin(symid), symax(symid);
  G4double szmin(szmid), szmax(szmid);
  G4Transform3D transform;
  switch (logoDirection) {
  case G4Scale::x:
    sxmin = sxmid - halfHeight;
    sxmax = sxmid + halfHeight;
    break;
  case G4Scale::y:
    symin = symid - halfHeight;
    symax = symid + halfHeight;
    transform = G4RotateZ3D(halfpi);
    break;
  case G4Scale::z:
    szmin = szmid - halfHeight;
    szmax = szmid + halfHeight;
    transform = G4RotateY3D(halfpi);
    break;
  }
  transform = G4Translate3D(sxmid,symid,szmid) * transform;

  G4VisAttributes visAtts(G4Colour(red, green, blue));
  visAtts.SetForceSolid(true);         // Always solid.

  G4Logo* logo = new G4Logo(height,visAtts);
  G4VModel* model =
    new G4CallbackModel<G4VisCommandSceneAddLogo::G4Logo>(logo);
  model->SetGlobalDescription("G4Logo");
  model->SetGlobalTag("G4Logo");
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
	     << ", ";
      switch (logoDirection) {
      case G4Scale::x:
	G4cout << 'x';
	break;
      case G4Scale::y:
	G4cout << 'y';
	break;
      case G4Scale::z:
	G4cout << 'z';
	break;
      }
      G4cout << "-direction, added to scene \"" << currentSceneName << "\"";
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
  fHeight(height),
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
  const G4double s = h;  // Half height of square subtractor
  const G4double y8 = s; // Choose y of subtractor for outer slope.
  const G4double x8 = ((-s * d - dx * (yt - y8)) / dy) + xt;
  G4double y9 = s; // Choose y of subtractor for inner slope.
  G4double x9 = ((-(s - w) * d - dx * (yt - y8)) / dy) + xt;
  // But to get inner, we make a triangle translated by...
  const G4double xtr = s - f1, ytr = -s - f2 -w;
  x9 += xtr; y9 += ytr;

  // G...
  G4Tubs tG("tG",ri,ro,d2,0.15*pi,1.85*pi);
  G4Box bG("bG",w2,ro2,d2);
  G4UnionSolid logoG("logoG",&tG,&bG,G4Translate3D(ri+w2,-ro2,0.));
  fpG = logoG.CreatePolyhedron();
  fpG->SetVisAttributes(&fVisAtts);
  fpG->Transform(G4Translate3D(-0.55*h,0.,0.));

  // 4...
  G4Box b1("b1",h2,h2,d2);
  G4Box bS("bS",s,s,d2+e);  // Subtractor.
  G4Box bS2("bS2",s,s,d2+2.*e);  // 2nd Subtractor.
  G4SubtractionSolid s1("s1",&b1,&bS,G4Translate3D(f1-s,f2-s,0.));
  G4SubtractionSolid s2("s2",&s1,&bS,G4Translate3D(f1+s+w,f2-s,0.));
  G4SubtractionSolid s3("s3",&s2,&bS,G4Translate3D(f1+s+w,f2+s+w,0.));
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4PSHitsModel* model = new G4PSHitsModel(newValue);
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
  fpCommand -> SetGuidance (G4Scale::GetGuidanceString());
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("length", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("direction", 's', omitable = true);
  parameter->SetGuidance ("'x', 'y' or 'z' - otherwise defaults to 'x'.");
  parameter->SetDefaultValue ("x");
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
  parameter =  new G4UIparameter ("auto|manual", 's', omitable = true);
  parameter->SetGuidance
    ("Automatic placement or manual placement at (xmid,ymid,zmid).");
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4double userLength, red, green, blue, xmid, ymid, zmid;
  G4String userLengthUnit, direction, auto_manual, positionUnit;
  std::istringstream is (newValue);
  is >> userLength >> userLengthUnit >> direction
     >> red >> green >> blue
     >> auto_manual
     >> xmid >> ymid >> zmid >> positionUnit;

  G4double length = userLength * G4UIcommand::ValueOf(userLengthUnit);
  G4double unit = G4UIcommand::ValueOf(positionUnit);
  xmid *= unit; ymid *= unit; zmid *= unit;

  std::ostringstream oss;
  oss << userLength << ' ' << userLengthUnit;
  G4String annotation(oss.str());

  G4Scale::Direction scaleDirection (G4Scale::x);
  if (direction(0) == 'y') scaleDirection = G4Scale::y;
  if (direction(0) == 'z') scaleDirection = G4Scale::z;

  G4bool autoPlacing = false; if (auto_manual == "auto") autoPlacing = true;
  // Parameters read and interpreted.

  // Useful constants, etc...
  const G4double halfLength(length / 2.);
  const G4double comfort(0.01);  // 0.15 seems too big.  0.05 might be better.
  const G4double onePlusComfort(1. + comfort);
  const G4double freeLengthFraction (1. + 2. * comfort);

  const G4VisExtent& sceneExtent = pScene->GetExtent();  // Existing extent.
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
    if (freeLengthFraction * (xmax - xmin) < length) room = false; break;
  case G4Scale::y:
    if (freeLengthFraction * (ymax - ymin) < length) room = false; break;
  case G4Scale::z:
    if (freeLengthFraction * (zmax - zmin) < length) room = false; break;
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
		false, xmid, ymid, zmid);
  G4VisAttributes* pVisAttr = new G4VisAttributes(G4Colour(red, green, blue));
  // Created of the heap because it needs a long lifetime.  This is a
  // mess.  The model determines the life but the vis atttributes are
  // associated with the scale.  There's no way of knowing when to
  // delete the vis atttributes!!!
  scale.SetVisAttributes(pVisAttr);
  G4VModel* model = new G4ScaleModel(scale);

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
  transform = G4Translate3D(sxmid,symid,szmid) * transform;
  //////////  G4VisExtent scaleExtent(sxmin, sxmax, symin, symax, szmin, szmax);


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
  fpCommand -> SetGuidance
    ("Adds text to current scene.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("x", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  parameter->SetGuidance ("x");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("y", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  parameter->SetGuidance ("y");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("z", 'd', omitable = true);
  parameter->SetDefaultValue (0);
  parameter->SetGuidance ("z");
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4bool smooth = false, rich = false;
  if (newValue.find("smooth") != std::string::npos) smooth = true;
  if (newValue.find("rich") != std::string::npos) rich = true;

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  G4PropagatorInField* propagatorInField =
    G4TransportationManager::GetTransportationManager()->
    GetPropagatorInField();
  propagatorInField->SetTrajectoryFilter(0); // Switch off smooth trajectories.
  static G4IdentityTrajectoryFilter auxiliaryPointsFilter;
  G4String defaultTrajectoryType;
  G4bool i_mode_found = false;
  G4int i_mode = 0;
  if (smooth && rich) {
    UImanager->ApplyCommand("/tracking/storeTrajectory 3");
    propagatorInField->SetTrajectoryFilter(&auxiliaryPointsFilter);
    defaultTrajectoryType = "G4RichTrajectory configured for smooth steps";
  } else if (smooth) {
    UImanager->ApplyCommand("/tracking/storeTrajectory 2");
    propagatorInField->SetTrajectoryFilter(&auxiliaryPointsFilter);
    defaultTrajectoryType = "G4SmoothTrajectory";
  } else if (rich) {
    UImanager->ApplyCommand("/tracking/storeTrajectory 3");
    defaultTrajectoryType = "G4RichTrajectory";
  } else {
    if (!newValue.empty()) {
      std::istringstream iss(newValue);
      iss >> i_mode;
      if (iss) {
	i_mode_found = true;
	if (verbosity >= G4VisManager::warnings) {
	  G4cout <<
  "WARNING: Integer parameter " << i_mode << " found."
  "\n  DEPRECATED - its use in this command will be removed at a future major"
  "\n  release.  Use \"/vis/modeling/trajectories\" commands."
		 << G4endl;
	}
      } else {
	if (verbosity >= G4VisManager::errors) {
	  G4cout << "ERROR: Unrecognised parameter \"" << newValue << "\""
	    "\n  No action taken."
		 << G4endl;
	}
	return;
      }
    }
    UImanager->ApplyCommand("/tracking/storeTrajectory 1");
    defaultTrajectoryType = "G4Trajectory";
  }
  UImanager->SetVerboseLevel(keepVerbose);

  if (rich) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"Attributes available for modeling and filtering with"
	"\n\"/vis/modeling/trajectories/create/drawByAttribute\" and"
	"\n\"/vis/filtering/trajectories/create/attributeFilter\" commands:\n"
	     << G4RichTrajectory().GetAttDefs()
	     << G4RichTrajectoryPoint().GetAttDefs();
    }
  }

  G4TrajectoriesModel* model = 0;
  if (i_mode_found) {
    model = new G4TrajectoriesModel(i_mode);
  } else {
    model = new G4TrajectoriesModel();
  }
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
  fpCommand = new G4UIcommand("/vis/scene/add/userAction",this);
  fpCommand -> SetGuidance
    ("Add Vis User Action, if any, to current scene.");
  fpCommand -> SetGuidance
    ("Optional arguments define the extent of the callback drawing.  You may"
     "\nnot need this if the extent has been defined in the original"
     "\nSetUserAction or is defined by other components of the scene.  But if"
     "\nthe user action is the only component of the scene, you will certainly"
     "\nneed to set the extent either in SetUserAction or here.  A scene must"
     "\nhave an extent one way or another so that the viewer can calculate"
     "\nhow to point the camera.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("xmin", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("xmax", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("ymin", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("ymax", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("zmin", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("zmax", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("cm");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSceneAddUserAction::~G4VisCommandSceneAddUserAction () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddUserAction::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddUserAction::SetNewValue (G4UIcommand*,
						    G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4bool warn = verbosity >= G4VisManager::warnings;

  G4VUserVisAction* visAction = fpVisManager->GetUserAction();
  if (!visAction) {
    if (warn) {
      G4cout <<	"WARNING: No User Vis Action registered." << G4endl;
    }
    return;
  }

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4String unitString;
  G4double xmin, xmax, ymin, ymax, zmin, zmax;
  std::istringstream is (newValue);
  is >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax >> unitString;
  G4double unit = G4UIcommand::ValueOf(unitString);
  xmin *= unit; xmax *= unit;
  ymin *= unit; ymax *= unit;
  zmin *= unit; zmax *= unit;
  G4VisExtent commandExtent(xmin,xmax,ymin,ymax,zmin,zmax);

  G4VisExtent extent;
  if (commandExtent.GetExtentRadius() > 0.) {
    extent = commandExtent;
  } else if (fpVisManager->GetUserActionExtent().GetExtentRadius() > 0.) {
    extent = fpVisManager->GetUserActionExtent();
  } else {
    if (warn) {
      G4cout <<	"WARNING: User Vis Action extent is null." << G4endl;
    }
  }

  G4VModel* model = new G4CallbackModel<G4VUserVisAction>(visAction);
  model->SetGlobalDescription("Vis User Action");
  model->SetGlobalTag("Vis User Action");
  model->SetExtent(extent);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful && verbosity >= G4VisManager::confirmations) {
    G4cout << "User Vis Action added to scene \""
	   << currentSceneName << "\"";
    if (verbosity >= G4VisManager::parameters) {
      G4cout << "\n  with extent " << extent;
    }
    G4cout << G4endl;
  }
  UpdateVisManagerScene (currentSceneName);
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
     "\n\"vis/scene/add/volume ! ! ! -box km 0 1 0 1 0 1\" will draw the world"
     "\nwith the positive octant cut away."); 
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
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
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

  G4VPhysicalVolume* world = *(transportationManager->GetWorldsIterator());

  if (!world) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<
	"ERROR: G4VisCommandSceneAddVolume::SetNewValue:"
	"\n  No world.  Maybe the geometry has not yet been defined."
	"\n  Try \"/run/initialize\""
	     << G4endl;
    }
    return;
  }

  const std::vector<G4VModel*>& rdModelList = pScene -> GetRunDurationModelList();
  std::vector<G4VModel*>::const_iterator i;
  for (i = rdModelList.begin(); i != rdModelList.end(); ++i) {
    if ((*i) -> GetGlobalDescription().find("G4PhysicalVolumeModel")
	!= std::string::npos) {
      if (((G4PhysicalVolumeModel*)(*i)) -> GetTopPhysicalVolume () == world) break;
    }
  }
  if (i != rdModelList.end()) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout << "WARNING: There is already a volume, \""
             << (*i) -> GetGlobalDescription()
             << "\",\n in the run-duration model list of scene \""
             << pScene -> GetName()
             << "\".\n To get a clean scene:"
	     << "\n  /vis/drawVolume " << name
	     << "\n or"
	     << "\n  /vis/scene/create"
	     << "\n  /vis/scene/add/volume " << name
	     << "\n  /vis/sceneHandler/attach"
	     << "\n (and also, if necessary, /vis/viewer/flush)"
             << G4endl;
    }
    return;
  }

  std::vector<G4PhysicalVolumeModel*> models;
  std::vector<G4VPhysicalVolume*> foundVolumes;
  G4VPhysicalVolume* foundWorld = 0;
  std::vector<G4int> foundDepths;
  std::vector<G4Transform3D> transformations;

  if (name == "world") {

    models.push_back
      (new G4PhysicalVolumeModel (world, requestedDepthOfDescent));
    foundVolumes.push_back(world);
    foundDepths.push_back(0);
    transformations.push_back(G4Transform3D());

  } else if (name == "worlds") {

    size_t nWorlds = transportationManager->GetNoWorlds();
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
    
    size_t nWorlds = transportationManager->GetNoWorlds();
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
	foundDepths.push_back(searchScene.GetFoundDepth());
	transformations.push_back(searchScene.GetFoundTransformation());
	break;
      }
    }

    if (foundVolumes.size()) {
      for (size_t i = 0; i < foundVolumes.size(); ++i) {
	models.push_back
	  (new G4PhysicalVolumeModel
	   (foundVolumes[i], requestedDepthOfDescent, transformations[i]));
      }
    } else {
      if (verbosity >= G4VisManager::errors) {
	G4cout << "ERROR: Volume \"" << name << "\"";
	if (copyNo >= 0) {
	  G4cout << ", copy no. " << copyNo << ",";
	}
	G4cout << " not found." << G4endl;
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
