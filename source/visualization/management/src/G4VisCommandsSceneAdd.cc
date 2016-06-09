//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VisCommandsSceneAdd.cc,v 1.59 2005/11/22 17:17:37 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
// /vis/scene commands - John Allison  9th August 1998

#include "G4VisCommandsSceneAdd.hh"

#include "G4VisManager.hh"
#include "G4TransportationManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4ModelingParameters.hh"
#include "G4HitsModel.hh"
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
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4Tokenizer.hh"
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
      G4bool successful(false);
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
  G4bool booleans, voxels, readout;
  std::istringstream is (newValue);
  is >> name >> requestedDepthOfDescent >>  booleans >> voxels >> readout;

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
	     << " not found in logical volume Store." << G4endl;
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
    return;
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
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
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
  G4bool warn(verbosity >= G4VisManager::warnings);

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

  G4bool autoPlacing (false); if (auto_manual(0) == 'a') autoPlacing = true;
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
  G4bool worried(false);
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
  G4bool room (true);
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

  G4Logo* logo = new G4Logo(height,visAtts,fpVisManager);
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
(G4double height, const G4VisAttributes& visAtts, G4VisManager* pVisManager):
  fHeight(height),
  fVisAtts(visAtts),
  fpVisManager(pVisManager)
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
  (const G4Transform3D& transform) {
  fpVisManager->Draw(*fpG,transform);
  fpVisManager->Draw(*fp4,transform);
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
  G4bool warn(verbosity >= G4VisManager::warnings);

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

  G4bool autoPlacing (false); if (auto_manual(0) == 'a') autoPlacing = true;
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
  G4bool worried(false);
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
  G4bool room (true);
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
  G4bool warn(verbosity >= G4VisManager::warnings);

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
  fpCommand = new G4UIcmdWithAnInteger
    ("/vis/scene/add/trajectories", this);
  fpCommand -> SetGuidance
    ("Adds trajectories to current scene.");
  fpCommand -> SetGuidance
    ("Causes trajectories, if any, to be drawn at the end of processiing an"
     "\nevent. The drawing mode is an integer that is passed to the"
     "\nDrawTrajectory method.  The default implementation in G4VTrajectory,"
     "\nif drawing-mode >= 0, draws the trajectory as a polyline (blue for"
     "\npositive, red for negative, green for neutral) and, if"
     "\ndrawing-mode != 0, draws markers of screen size"
     "\nstd::abs(drawing-mode)/1000 pixels at each step and auxiliary point,"
     "\nif any.  So drawing-mode = 5000 is a good choice."
     "\nEnable storing with \"/tracking/storeTrajectory 1\"."
     "\nSee also \"/vis/scene/endOfEventAction\".");
  fpCommand -> SetParameterName ("drawing-mode", omitable = true);
  fpCommand -> SetDefaultValue (0);
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
  G4bool warn(verbosity >= G4VisManager::warnings);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<	"ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4int drawingMode;
  std::istringstream is (newValue);
  is >> drawingMode;
  G4TrajectoriesModel* model = new G4TrajectoriesModel(drawingMode);
  const G4String& currentSceneName = pScene -> GetName ();
  pScene -> AddEndOfEventModel (model, warn);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose = 0;
  if (keepVerbose >= 2 ||
      fpVisManager->GetVerbosity() >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  UImanager->ApplyCommand("/tracking/storeTrajectory 1");
  UImanager->SetVerboseLevel(keepVerbose);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Trajectories will be drawn with mode "
	   << drawingMode
	   << " in scene \""
	   << currentSceneName << "\"."
	   << G4endl;
  }
  if (verbosity >= G4VisManager::warnings) {
    G4cout << "WARNING: \"/tracking/storeTrajectory 1\" has been executed."
	   << G4endl;
  }
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
  G4bool warn(verbosity >= G4VisManager::warnings);

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
  pScene -> AddRunDurationModel (model, warn);
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "User Vis Action added to scene \""
	   << currentSceneName << "\"";
    if (verbosity >= G4VisManager::parameters) {
      G4cout << "\n  with extent " << extent;
    }
    G4cout << G4endl;
  }
}

////////////// /vis/scene/add/volume ///////////////////////////////////////

G4VisCommandSceneAddVolume::G4VisCommandSceneAddVolume () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/volume", this);
  fpCommand -> SetGuidance 
   ("Adds a physical volume to current scene, with optional clipping volume.");
  fpCommand -> SetGuidance 
    ("If physical-volume-name is \"world\" (the default), the top of the"
     "\ntracking tree is used (GetNavigatorForTracking()->GetWorldVolume())."
     "\nOtherwise a search of the tracking tree is made, taking the first"
     "\nmatching occurence only.  To see a representation of the geometry"
     "\nhierarchy of the tracking tree, try \"/vis/drawTree\" or one of the"
     "\ndriver/browser combinations that have the required functionality,"
     "\ne.g., HepRepFile/XML with the WIRED3/4 browser.");
  fpCommand -> SetGuidance
    ("If clip-volume-type is specified, the subsequent parameters are used to"
     "\nto define a clipping volume.  For example,"
     "\n\"vis/scene/add/volume ! ! ! box km 0 1 0 1 0 1\" will draw the world"
     "\nwith the positive octant cut away."); 
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
  parameter -> SetDefaultValue ("none");
  parameter -> SetGuidance
    ("For \"box\", the parameters are xmin,xmax,ymin,ymax,zmin,zmax."
     // "\n Only \"box\" is programmed at present."); No '\n' for GAG (temp).
     " Only \"box\" is programmed at present.");
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
  G4bool warn(verbosity >= G4VisManager::warnings);

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
  G4double unit = G4UIcommand::ValueOf(parameterUnit);
  param1 *= unit; param2 *= unit; param3 *= unit;
  param4 *= unit; param5 *= unit; param6 *= unit;

  G4VPhysicalVolume* world =
    G4TransportationManager::GetTransportationManager ()
    -> GetNavigatorForTracking () -> GetWorldVolume ();
  G4PhysicalVolumeModel* model = 0;
  G4VPhysicalVolume* foundVolume = 0;
  G4int foundDepth = 0;

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
  }

  if (name == "world") {
    if (world) {
      model = new G4PhysicalVolumeModel (world,
					 requestedDepthOfDescent);
      foundVolume = world;
    }
    else {
      if (verbosity >= G4VisManager::errors) {
	G4cout << "ERROR: G4VisCommandSceneAddVolume::SetNewValue:"
	       << "\n  No world - shouldn't happen if G4ApplicationState is"
	       << " being properly noted!!" << G4endl;
      }
      return;
    }
  }
  else {
 
    // Create search scene, model and modeling parameters with
    // long-enough life...
    G4PhysicalVolumeSearchScene searchScene (name, copyNo);
    G4PhysicalVolumeModel searchModel (world);  // Default - unlimited depth.
    G4ModelingParameters mp;  // Default - no culling.
    searchModel.SetModelingParameters (&mp);

    // Initiate search...
    searchModel.DescribeYourselfTo (searchScene);

    // OK, what have we got...?
    foundVolume = searchScene.GetFoundVolume ();
    foundDepth = searchScene.GetFoundDepth ();
    const G4Transform3D&
      transformation = searchScene.GetFoundTransformation ();

    if (foundVolume) {
      model = new G4PhysicalVolumeModel (foundVolume,
					 requestedDepthOfDescent,
					 transformation);
    }
    else {
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
    G4Box clippingBox("_clipping_box",dX,dY,dZ);
    G4Polyhedron* clippingPolyhedron = new G4PolyhedronBox(dX,dY,dZ);
    // Created on the heap and left there.
    clippingPolyhedron->Transform(G4Translate3D(x0,y0,z0));
    model->SetClippingPolyhedron(clippingPolyhedron);
  }  // If any other shape consider NumberOfRotationSides!!!!!!!!!!!

  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "First occurrence of \""
	     << foundVolume -> GetName ()
	     << "\"";
      if (copyNo >= 0) {
	G4cout << ", copy no. " << copyNo << ",";
      }
      G4cout << " found at depth " << foundDepth
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
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
  UpdateVisManagerScene (currentSceneName);
}
