// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManMessDraw.cc,v 1.7 2001-02-23 15:43:29 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager Messenger - John Allison 22nd July 1996.
// Draw sub-menu.

#include "G4VisManMessenger.hh"

#include "G4VisManager.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "G4VGlobalFastSimulationManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4FlavoredParallelWorldModel.hh"
#include "G4ModelingParameters.hh"


void G4VisManMessenger::AddCommandDraw () {

  G4UIcommand* command;
  G4UIparameter* param;

  ////////////////////////////////////  /vis~/draw/axes  /////
  //draw \hline
  //draw /vis~/draw/axes & $x_0$, $y_0$, $z_0$, length, unit &
  //draw Draws axes at $(x_0, y_0, z_0)$ of given length. \\%
  command = new G4UIcommand ("/vis~/draw/axes", this);
  command -> SetGuidance
    (
     "Draws axes at (x0, y0, z0) of given length."
     );
  param   =  new G4UIparameter ("x0", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("y0", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("z0", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("length", 'd', true);
  param   -> SetDefaultValue (10000.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("unit", 's', true);
  param   -> SetDefaultValue  ("mm");
  param   -> SetGuidance      ("mm, cm, or m.");
  command -> SetParameter     (param);
  fCommandList.push_back (command);

  ////////////////////////////////////  /vis~/draw/text  /////
  //draw \hline
  //draw /vis~/draw/text & $x$, $y$, $z$, unit,
  //draw font\_size, x\_offset, y\_offset, text &
  //draw Draws text at $(x, y, z)$ with given parameters.
  //draw Font size and offsets in pixels. \\%
  command = new G4UIcommand ("/vis~/draw/text", this);
  command -> SetGuidance 
    ("Draws at (x, y, z) unit font_size x_offset y_offset text.");
  command -> SetGuidance
    ("Font size and offsets in pixels.");
  param   =  new G4UIparameter ("x", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("y", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("z", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("unit", 's', true);
  param   -> SetDefaultValue  ("mm");
  param   -> SetGuidance      ("mm, cm, or m.");
  command -> SetParameter     (param);
  param   =  new G4UIparameter ("font_size", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("x_offset", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("y_offset", 'd', true);
  param   -> SetDefaultValue (0.);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("text", 's', true);
  param   -> SetDefaultValue ("text");
  command -> SetParameter (param);
  fCommandList.push_back (command);

  //
  ////////////////////////////////////  /vis~/draw/ghosts  /////
  command = new G4UIcommand ("/vis~/draw/ghosts", this);
  command -> SetGuidance
    (
     "Adds to the current scene all ghost volumes"
     );
  param   =  new G4UIparameter ("ParticleName", 's', true);
  param   -> SetDefaultValue ("all");
  command -> SetParameter (param);
  fCommandList.push_back (command);


  /********************************  UNDER DEVELOPMENT ?????????????
  //////////////////////////  /vis~/draw/volume  /////
  //??draw \hline
  //??draw /vis~/draw/volume & name &
  //??draw Clears the current scene, makes a new one consisting of the
  //??draw physical volume ``name'' and
  //??draw draws it in the current view, marking all other views of this scene
  //??draw as needing refreshing.  \\%
  command = new G4UIcommand ("/vis~/draw/volume", this);
  command -> SetGuidance
    (
     "Clears the current scene, makes a new one consisting of the "
     "physical volume \"name\" and "
     "draws it in the current view, marking all other views of this scene "
     "as needing refreshing."
     );
  param   =  new G4UIparameter ("Physical volume name", 's', true);
  param   -> SetDefaultValue  ("__none__");
  command -> SetParameter     (param);
  fCommandList.push_back (command);
  **************************************************************/

}

void G4VisManMessenger::DoCommandDraw (const G4String& commandPath,
				       G4String& newValues) {

  /////////////////////////////////////////  /vis~/draw/axes  ////
  if (commandPath == "/vis~/draw/axes") {
    if (ViewValid ()) {
      G4double x0, y0, z0;
      G4double length ;
      G4String unitString;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString);
      is >> x0 >> y0 >> z0 >> length >> unitString;

      G4double unit = G4UnitDefinition::GetValueOf (unitString);

      x0 *= unit; y0 *= unit; z0 *= unit; length *= unit;

      G4Polyline x_axis, y_axis, z_axis;

      G4Colour cx(1.0, 0.0, 0.0); // color for x-axis
      G4Colour cy(0.0, 1.0, 0.0); // color for y-axis
      G4Colour cz(0.0, 0.0, 1.0); // color for z-axis

      G4VisAttributes ax(cx);     // VA for x-axis
      G4VisAttributes ay(cy);     // VA for y-axis
      G4VisAttributes az(cz);     // VA for z-axis

		//----- Draw x-axis
      x_axis.SetVisAttributes(&ax);
      x_axis.append (G4Point3D (x0*unit,y0,z0));
      x_axis.append ( G4Point3D ( (x0 +length) , y0, z0 ) );
      fpVMan -> Draw (x_axis);
      x_axis.clear ();

		//----- Draw y-axis
      y_axis.SetVisAttributes(&ay);
      y_axis.append (G4Point3D (x0,y0,z0));
      y_axis.append ( G4Point3D ( x0 , (y0 +length ) , z0 ) );
      fpVMan -> Draw (y_axis);
      y_axis.clear ();

		//----- Draw z-axis
      z_axis.SetVisAttributes(&az);
      z_axis.append (G4Point3D (x0,y0,z0));
      z_axis.append ( G4Point3D ( x0 , y0 , ( z0 +length ) ) );
      fpVMan -> Draw (z_axis);
      z_axis.clear ();

    }
  }

  /////////////////////////////////////////  /vis~/draw/text  ////
  if (commandPath == "/vis~/draw/text") {
    if (ViewValid ()) {
      G4double x, y, z;
      G4double font_size, x_offset, y_offset;
      G4String unitString, textString;
      const char* aString = newValues;
      G4std::istrstream is((char*) aString);
      is >> x >> y >> z >> unitString
	 >> font_size >> x_offset >> y_offset >> textString;
      G4double unit = G4UnitDefinition::GetValueOf (unitString);
      x *= unit; y *= unit; z *= unit;
      G4Text g4Text (textString, G4Point3D (x, y, z));
      g4Text.SetScreenSize (font_size);
      g4Text.SetOffset (x_offset, y_offset);
      fpVMan -> Draw (g4Text);
    }
  }

  //////////////////////////  /vis~/draw/volume  /////
  if (commandPath == "/vis~/draw/volume") {
    G4VisManager::PrintCommandDeprecation
      ("Use \"/vis/scene/add/volume\" or \"/vis/drawVolume\".");
     // Find physical_volume.
    G4String name;
    const char* aString = newValues;
    G4std::istrstream is((char*) aString) ; is >> name;

    G4cout << "/vis~/draw/volume: requested name is " << name;
    G4cout << "\n...but this command is \"UNDER DEVELOPMENT\"!";
    G4cout << G4endl;
  }
  //////////////////////////  /vis~/draw/ghosts  /////
  if (commandPath == "/vis~/draw/ghosts") {
    G4VisManager::PrintCommandDeprecation
      ("Use \"/vis/scene/add/ghosts\".");
    // preliminaries
    G4VisManager* theVisManager;
    if(!(theVisManager= G4VisManager::GetInstance())) {
      G4cout<< "grafic system not yet avaliable." << G4endl;
      return;
    }
    G4VGlobalFastSimulationManager* theGlobalFastSimulationManager;
    if(!(theGlobalFastSimulationManager = 
	 G4VGlobalFastSimulationManager::GetConcreteInstance ())){
      G4cout<< "WARNING: none G4GlobalFastSimulationManager" << G4endl;
      return;
    }
    G4ParticleTable* theParticleTable=G4ParticleTable::GetParticleTable();
    
    G4Scene* currentScene = theVisManager -> GetCurrentScene ();
    // Ok, lets look on the parameters
    if(newValues=="all") {
      G4VFlavoredParallelWorld* CurrentFlavoredWorld;
      for (G4int iParticle=0; iParticle<theParticleTable->entries(); 
	   iParticle++)
	if(CurrentFlavoredWorld=theGlobalFastSimulationManager->
	   GetFlavoredWorldForThis(theParticleTable->
				   GetParticle(iParticle)))
	  currentScene -> AddRunDurationModel
	    (new G4FlavoredParallelWorldModel (CurrentFlavoredWorld));
      G4cout << "Ghosts added to the Scene, refresh the view to see it."
	     << G4endl;
      return;
    }
    G4ParticleDefinition* currentParticle = 
      theParticleTable->FindParticle(newValues);
    if (currentParticle == NULL) {
      G4cout << newValues << ": not found this particle name!" << G4endl;
      return;
    }
    G4VFlavoredParallelWorld* worldForThis;
    if(worldForThis=theGlobalFastSimulationManager->
       GetFlavoredWorldForThis(currentParticle)) {
      currentScene -> AddRunDurationModel
	(new G4FlavoredParallelWorldModel (worldForThis));
      G4cout << "Ghosts added to the Scene, refresh the view to see it."
	     << G4endl;
    }
    else G4cout << "There are no ghosts for "<<newValues<<G4endl;
  }
}
