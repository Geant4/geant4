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
// $Id: G4VisCommandsSceneAdd.cc,v 1.26 2001-08-14 18:32:45 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4ios.hh"
#include "g4std/strstream"

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
  fpCommand -> SetGuidance
    ("/vis/scene/add/axes [<x0>] [<y0>] [<z0>] [<length>] [<unit>]");
    ("Default: 0 0 0 1 m");
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
  parameter->SetGuidance      ("mm, cm, or m.");
  fpCommand->SetParameter     (parameter);
}

G4VisCommandSceneAddAxes::~G4VisCommandSceneAddAxes () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddAxes::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddAxes::SetNewValue (G4UIcommand* command,
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

  G4String unitString;
  G4double x0, y0, z0, length;
  const char* s = newValue;
  G4std::istrstream is ((char*)s);
  is >> x0 >> y0 >> z0 >> length >> unitString;

  G4double unit = ValueOf(unitString);
  x0 *= unit; y0 *= unit; z0 *= unit; length *= unit;

  G4VModel* model = new G4AxesModel(x0, y0, z0, length);
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
    ("/vis/scene/add/ghosts [<particle-name>]");
  fpCommand -> SetGuidance
    ("Adds ghost volumes (G4FlavoredParallelWorld) to the current scene.");
  fpCommand -> SetGuidance
    ("Selects by particle (default = \"all\").");
  fpCommand -> SetParameterName ("particle", omitable = true);
  fpCommand -> SetDefaultValue ("all");
}

G4VisCommandSceneAddGhosts::~G4VisCommandSceneAddGhosts () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddGhosts::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneAddGhosts::SetNewValue (G4UIcommand* command,
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
  const G4String& currentSceneName = pScene -> GetName ();

  G4VGlobalFastSimulationManager* theGlobalFastSimulationManager;
  if(!(theGlobalFastSimulationManager = 
       G4VGlobalFastSimulationManager::GetConcreteInstance ())){
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: no G4GlobalFastSimulationManager" << G4endl;
    }
    return;
  }
  
  G4ParticleTable* theParticleTable=G4ParticleTable::GetParticleTable();
  
  if(newValue=="all") {
    G4VFlavoredParallelWorld* CurrentFlavoredWorld;
    G4bool successful;
    for (G4int iParticle=0; iParticle<theParticleTable->entries(); 
	 iParticle++)
      CurrentFlavoredWorld = theGlobalFastSimulationManager->
	GetFlavoredWorldForThis(theParticleTable->GetParticle(iParticle));
    if(CurrentFlavoredWorld)
      successful = pScene -> AddRunDurationModel
	(new G4FlavoredParallelWorldModel (CurrentFlavoredWorld), warn);
    if (successful) {
      if (verbosity >= G4VisManager::confirmations) {
	G4cout << "Ghosts have been added to scene \""
	       << currentSceneName << "\"."
	       << G4endl;
      }
    }
    else {
      G4VisCommandsSceneAddUnsuccessful(verbosity);
      return;
    }
  }
  
  G4ParticleDefinition* currentParticle = 
    theParticleTable->FindParticle(newValue);
  if (currentParticle == NULL) {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: \"" << newValue
	     << "\": not found this particle name!" << G4endl;
    }
    return;
  }

  G4VFlavoredParallelWorld* worldForThis =
    theGlobalFastSimulationManager->GetFlavoredWorldForThis(currentParticle);
  if(worldForThis) {
    G4bool successful = pScene -> AddRunDurationModel
      (new G4FlavoredParallelWorldModel (worldForThis), warn);
    if (successful) {
      if (verbosity >= G4VisManager::confirmations) {
	G4cout << "Ghosts have been added to scene \""
	       << currentSceneName << "\"."
	       << G4endl;
      }
    }
  }
  else {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: There are no ghosts for \""<<newValue<<"\""<<G4endl;
      G4VisCommandsSceneAddUnsuccessful(verbosity);
    }
  }
  UpdateVisManagerScene (currentSceneName);
}


////////////// /vis/scene/add/hits ///////////////////////////////////////

G4VisCommandSceneAddHits::G4VisCommandSceneAddHits () {
  fpCommand = new G4UIcmdWithoutParameter ("/vis/scene/add/hits", this);
  fpCommand -> SetGuidance
    ("Adds hits to current scene.");
  fpCommand -> SetGuidance
    ("Hits are drawn at end of event when the scene in which"
     " they are added is current.");
}

G4VisCommandSceneAddHits::~G4VisCommandSceneAddHits () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddHits::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneAddHits::SetNewValue (G4UIcommand* command,
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
  fpCommand -> SetGuidance
    ("/vis/scene/add/logicalVolume <logical-volume-name>"
     " [<depth-of-descending>]");
  fpCommand -> SetGuidance ("Adds a logical volume to the current scene.");
  fpCommand -> SetGuidance
    ("1st parameter: volume name.");
  //  fpCommand -> SetGuidance  // Not implemented - should be in geom?
  //    ("               \"list\" to list all volumes.");
  fpCommand -> SetGuidance
    ("2nd parameter: depth of descending geometry hierarchy (default 1).");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("volume", 's', omitable = false);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("depth", 'i', omitable = true);
  parameter -> SetDefaultValue (1);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddLogicalVolume::~G4VisCommandSceneAddLogicalVolume () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddLogicalVolume::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddLogicalVolume::SetNewValue (G4UIcommand* command,
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
  const char* s = newValue;
  G4std::istrstream is ((char*)s);
  is >> name >> requestedDepthOfDescent;

  G4LogicalVolumeStore *pLVStore = G4LogicalVolumeStore::GetInstance();
  int nLV = pLVStore -> size ();
  int iLV;
  G4LogicalVolume* pLV;
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

  G4VModel* model = new G4LogicalVolumeModel (pLV, requestedDepthOfDescent);
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Logical volume \"" << pLV -> GetName ()
	     << " with requested depth of descent "
	     << requestedDepthOfDescent
	     << ",\n  has been added to scene \"" << currentSceneName << "\"."
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
    ("/vis/scene/add/scale [<length> <length-unit>] [x|y|z] [<red>] [<green>] [<blue>] [auto|manual] [<xmid> <ymid> <zmid> <unit>]");
  fpCommand -> SetGuidance 
    ("Defaults: 1 m x 1 1 1 auto 0 0 0 m");
  fpCommand -> SetGuidance 
    ("Adds an annotated scale line to the current scene.");
  fpCommand -> SetGuidance 
    ("See G4Scale.hh for further description.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("length", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue ("m");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("direction", 's', omitable = true);
  parameter->SetDefaultValue ("x");
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("red", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("auto|manual", 's', omitable = true);
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

void G4VisCommandSceneAddScale::SetNewValue (G4UIcommand* command,
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

  G4double userLength, red, green, blue, xmid, ymid, zmid;
  G4String userLengthUnit, direction, auto_manual, positionUnit;
  G4std::istrstream is (newValue);
  is >> userLength >> userLengthUnit >> direction
     >> red >> green >> blue
     >> auto_manual
     >> xmid >> ymid >> zmid >> positionUnit;

  G4double length = userLength * ValueOf(userLengthUnit);
  G4double unit = ValueOf(positionUnit);
  xmid *= unit; ymid *= unit; zmid *= unit;

  char tempcharstring [50];
  G4std::ostrstream ost (tempcharstring, 50);
  ost << userLength << ' ' << userLengthUnit << G4std::ends;
  G4String annotation(tempcharstring);

  G4Scale::Direction scaleDirection (G4Scale::x);
  if (direction(0) == 'y') scaleDirection = G4Scale::y;
  if (direction(0) == 'z') scaleDirection = G4Scale::z;

  G4bool autoPlacing (false); if (auto_manual(0) == 'a') autoPlacing = true;
  // Parameters read and interpreted.

  // Useful constants, etc...
  const G4double halfLength(length / 2.);
  const G4double comfort(0.15);
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

  // Let's go ahead a construct a scale and a scale model...
  G4Scale scale(length, annotation, scaleDirection,
		autoPlacing, xmid, ymid, xmid);
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
  switch (scaleDirection) {
  case G4Scale::x:
    sxmin = sxmid - halfLength;
    sxmax = sxmid + halfLength;
    break;
  case G4Scale::y:
    symin = symid - halfLength;
    symax = symid + halfLength;
    break;
  case G4Scale::z:
    szmin = szmid - halfLength;
    szmax = szmid + halfLength;
    break;
  }
  G4VisExtent scaleExtent(sxmin, sxmax, symin, symax, szmin, szmax);
  model->SetExtent(scaleExtent);
  // This extent gets "added" to existing scene extent in
  // AddRunDurationModel below.

  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddRunDurationModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Scale of " << annotation
	     << " has been added to scene \"" << currentSceneName << "\"."
	     << G4endl;
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
    ("Adds text at (x, y, z) unit font_size x_offset y_offset text.");
  fpCommand -> SetGuidance
    ("Font size and offsets in pixels.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("x", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("y", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("z", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("unit", 's', omitable = true);
  parameter->SetDefaultValue  ("mm");
  parameter->SetGuidance      ("mm, cm, or m.");
  fpCommand->SetParameter     (parameter);
  parameter =  new G4UIparameter ("font_size", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("x_offset", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("y_offset", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter =  new G4UIparameter ("text", 's', omitable = true);
  parameter->SetDefaultValue ("text");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSceneAddText::~G4VisCommandSceneAddText () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddText::GetCurrentValue (G4UIcommand*) {
  return "";
}

void G4VisCommandSceneAddText::SetNewValue (G4UIcommand* command,
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

  G4String text, unitString;
  G4double x, y, z, font_size, x_offset, y_offset;
  const char* s = newValue;
  G4std::istrstream is ((char*)s);
  is >> x >> y >> z >> unitString >> font_size >> x_offset >> y_offset >> text;

  G4double unit = ValueOf(unitString);
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
  fpCommand = new G4UIcmdWithoutParameter
    ("/vis/scene/add/trajectories", this);
  fpCommand -> SetGuidance
    ("Adds trajectories to current scene.");
  fpCommand -> SetGuidance
    ("Trajectories are drawn at end of event when the scene in which"
     " they are added is current.");
}

G4VisCommandSceneAddTrajectories::~G4VisCommandSceneAddTrajectories () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddTrajectories::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneAddTrajectories::SetNewValue (G4UIcommand* command,
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

  G4TrajectoriesModel* model = new G4TrajectoriesModel;
  const G4String& currentSceneName = pScene -> GetName ();
  G4bool successful = pScene -> AddEndOfEventModel (model, warn);
  if (successful) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Trajectories will be drawn in scene \""
	     << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
}

////////////// /vis/scene/add/volume ///////////////////////////////////////

G4VisCommandSceneAddVolume::G4VisCommandSceneAddVolume () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/volume", this);
  fpCommand -> SetGuidance
    ("/vis/scene/add/volume [<physical-volume-name>] [<copy-no>] [<depth-of-descending>]");
  fpCommand -> SetGuidance ("Adds a physical volume to the current scene.");
  fpCommand -> SetGuidance ("Note: adds first occurence only.");
  fpCommand -> SetGuidance
    ("1st parameter: volume name (default \"world\").");
  //  fpCommand -> SetGuidance  // Not implemented - should be in geom?
  //    ("               \"list\" to list all volumes.");
  fpCommand -> SetGuidance
    ("2nd parameter: copy number (default -1)."
     "\n  If negative, first occurrence of physical-volume-name is selected.");
  fpCommand -> SetGuidance
    ("3rd parameter: depth of descending geometry hierarchy"
     " (default G4Scene::UNLIMITED (-1)).");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("volume", 's', omitable = true);
  parameter -> SetDefaultValue ("world");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("copy-no", 'i', omitable = true);
  parameter -> SetDefaultValue (-1);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("depth", 'i', omitable = true);
  parameter -> SetDefaultValue (G4Scene::UNLIMITED);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddVolume::~G4VisCommandSceneAddVolume () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddVolume::GetCurrentValue (G4UIcommand* command) {
  return "world 0 -1";
}

void G4VisCommandSceneAddVolume::SetNewValue (G4UIcommand* command,
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
  G4int copyNo;
  G4int requestedDepthOfDescent;
  const char* s = newValue;
  G4std::istrstream is ((char*)s);
  is >> name >> copyNo >> requestedDepthOfDescent;
  G4VPhysicalVolume* world =
    G4TransportationManager::GetTransportationManager ()
    -> GetNavigatorForTracking () -> GetWorldVolume ();
  G4PhysicalVolumeModel* model = 0;
  G4VPhysicalVolume* foundVolume = 0;
  G4int foundDepth = 0;

  if (name == "world") {
    if (world) {
      model = new G4PhysicalVolumeModel (world);
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
    G4PhysicalVolumeModel searchModel (world);
    G4ModelingParameters mp;
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
	     << ",\n  with ";
      if (requestedDepthOfDescent < 0) {
	G4cout << "unlimited (-1)";
      }
      else {
	G4cout << requestedDepthOfDescent;
      }
      G4cout << " further requested depth of descent"
	     << ",\n  has been added to scene \"" << currentSceneName << "\"."
	     << G4endl;
    }
  }
  else G4VisCommandsSceneAddUnsuccessful(verbosity);
  UpdateVisManagerScene (currentSceneName);
}
