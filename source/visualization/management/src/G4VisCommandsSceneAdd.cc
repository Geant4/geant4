// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsSceneAdd.cc,v 1.1 1999-01-07 16:15:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/scene commands - John Allison  9th August 1998

#include "G4VisCommandsSceneAdd.hh"

#include "G4VisManager.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ModelingParameters.hh"
#include "G4PhysicalVolumeSearchScene.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4FlavoredParallelWorldModel.hh"
#include "G4ApplicationState.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

////////////// /vis/scene/add/volume ///////////////////////////////////////

G4VisCommandSceneAddVolume::G4VisCommandSceneAddVolume () {
  G4bool omitable;
  fpCommand = new G4UIcommand ("/vis/scene/add/volume", this);
  fpCommand -> AvailableForStates (Idle, GeomClosed);
  fpCommand -> SetGuidance
    ("/vis/scene/add/volume [<physical-volume-name>] [<copy-no>] [<depth-of-descending>]");
  fpCommand -> SetGuidance ("Adds a physical volume to the current scene.");
  fpCommand -> SetGuidance
    ("1st parameter: volume name (default \"world\").");
  fpCommand -> SetGuidance
    ("2nd parameter: copy number (default 0).");
  fpCommand -> SetGuidance
    ("3rd parameter: depth of descending geometry hierarchy"
     " (default G4SceneData::UNLIMITED (-1)).");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("volume", 's', omitable = true);
  parameter -> SetDefaultValue ("world");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("copy-no", 'i', omitable = true);
  parameter -> SetDefaultValue (0);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter ("depth", 'i', omitable = true);
  parameter -> SetDefaultValue (G4SceneData::UNLIMITED);
  fpCommand -> SetParameter (parameter);
}

G4VisCommandSceneAddVolume::~G4VisCommandSceneAddVolume () {
  delete fpCommand;
}

G4String G4VisCommandSceneAddVolume::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneAddVolume::SetNewValue (G4UIcommand* command,
					      G4String newValue) {
  G4SceneDataObjectList& list = fpVisManager -> SetSceneDataObjectList ();
  if (list.isEmpty ()) {
    G4cout << "No scenes - please create one before adding anything."
	   << endl;
    return;
  }

  G4String name;
  G4int copyNo;
  G4int requestedDepthOfDescent;
  const char* s = newValue;
  istrstream is ((char*)s);
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
      G4cerr << "G4VisCommandSceneAddVolume::SetNewValue: *** ERROR ***"
	     << "\n  No world - shouldn't happen if G4ApplicationState is"
	     << " being properly handled!!" << endl;
    }
  }
  else {
    G4PhysicalVolumeSearchScene searchScene (name, copyNo);
    G4PhysicalVolumeModel searchModel (world);
    searchModel.DescribeYourselfTo (searchScene);
    foundVolume = searchScene.GetFoundVolume ();
    const G4Transform3D&
      transformation = searchScene.GetFoundTransformation ();
    foundDepth = searchScene.GetFoundDepth ();
    if (foundVolume) {
      model = new G4PhysicalVolumeModel (foundVolume,
					 requestedDepthOfDescent,
					 transformation);
    }
    else {
      G4cout << "Volume \"" << name << "\", copy no. " << copyNo
	     << " not found." << endl;
    }
  }

  if (model) {
    const G4String& currentSceneName =
      fpVisManager -> GetCurrentSceneData ().GetName ();
    (list [currentSceneName]).AddRunDurationModel (model);
    UpdateVisManagerSceneDataAndViewParameters (currentSceneName);
    G4cout << "First occurrence of \"" << foundVolume -> GetName ()
	   << "\", copy no. " << copyNo
	   << ", found at depth " << foundDepth
	   << ",\n  with further requested depth of descent "
	   << requestedDepthOfDescent
	   << ", has been added to scene \"" << currentSceneName << "\""
	   << endl;
  }
}

////////////// /vis/scene/add/ghosts ///////////////////////////////////////

G4VisCommandSceneAddGhosts::G4VisCommandSceneAddGhosts () {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString ("/vis/scene/add/ghosts", this);
  fpCommand -> AvailableForStates (Idle, GeomClosed);
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
  const G4String& currentSceneName =
    fpVisManager -> GetCurrentSceneData ().GetName ();

  G4SceneDataObjectList& list = fpVisManager -> SetSceneDataObjectList ();
  if (list.isEmpty ()) {
    G4cout << "No scenes - please create one before adding anything."
	   << endl;
    return;
  }

  G4GlobalFastSimulationManager* theGlobalFastSimulationManager;
  if(!(theGlobalFastSimulationManager = 
       G4GlobalFastSimulationManager::GetGlobalFastSimulationManager())){
    G4cout<< "WARNING: no G4GlobalFastSimulationManager" << endl;
    return;
  }

  G4ParticleTable* theParticleTable=G4ParticleTable::GetParticleTable();

  if(newValue=="all") {
    G4FlavoredParallelWorld* CurrentFlavoredWorld;
    for (G4int iParticle=0; iParticle<theParticleTable->entries(); 
	 iParticle++)
      if(CurrentFlavoredWorld=theGlobalFastSimulationManager->
	 GetFlavoredWorldForThis(theParticleTable->
				 GetParticle(iParticle)))
	(list [currentSceneName]).AddRunDurationModel
	  (new G4FlavoredParallelWorldModel (CurrentFlavoredWorld));
    UpdateVisManagerSceneDataAndViewParameters ();
    G4cout << "Ghosts added to the Scene, refresh the view to see it."
	   << endl;
    return;
  }

  G4ParticleDefinition* currentParticle = 
    theParticleTable->FindParticle(newValue);
  if (currentParticle == NULL) {
    G4cout << "\"" << newValue << "\": not found this particle name!" << endl;
    return;
  }

  G4FlavoredParallelWorld* worldForThis;
  if(worldForThis=theGlobalFastSimulationManager->
     GetFlavoredWorldForThis(currentParticle)) {
    (list [currentSceneName]).AddRunDurationModel
      (new G4FlavoredParallelWorldModel (worldForThis));
    UpdateVisManagerSceneDataAndViewParameters (currentSceneName);
    G4cout << "Ghosts added to the Scene, refresh the view to see it."
           << endl;
  }
  else G4cout << "There are no ghosts for \""<<newValue<<"\""<<endl;
}
