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

#include "G4OpenGLViewerMessenger.hh"

#include "G4OpenGLViewer.hh"
#include "G4OpenGLStoredViewer.hh"
#include "G4OpenGLStoredSceneHandler.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4VisManager.hh"
#include <sstream>

G4OpenGLViewerMessenger*
G4OpenGLViewerMessenger::fpInstance = 0;

G4OpenGLViewerMessenger*
G4OpenGLViewerMessenger::GetInstance()
{
  if (!fpInstance) fpInstance = new G4OpenGLViewerMessenger;
  return fpInstance;
}

G4OpenGLViewerMessenger::G4OpenGLViewerMessenger()
{
  G4bool omitable;

  fpDirectory = new G4UIdirectory("/vis/ogl/");
  fpDirectory->SetGuidance("G4OpenGLViewer commands.");

  fpCommandExport =
  new G4UIcommand("/vis/ogl/export", this);
  fpCommandExport->SetGuidance ("Export a screenshot of current OpenGL viewer");
  fpCommandExport->SetGuidance
  ("If name is \"\", filename and extension will have the current value");
  fpCommandExport->SetGuidance
  ("If name is \"toto.png\", set the name to \"toto\" and the format to \"png\".");
  fpCommandExport->SetGuidance
  ("If name is \"toto\", set the name to \"toto\" and the format to current format.");
  fpCommandExport->SetGuidance
  ("Will also add an incremented suffix at the end of the name, except if name is"
   "\nthe same as previous it will not reset the incremented suffix.");
  fpCommandExport->SetGuidance
  ("Setting size is available only on eps/pdf/svg/ps formats.");
  G4UIparameter* parameterExport;
  parameterExport = new G4UIparameter ("name", 's', omitable = true);
  parameterExport->SetDefaultValue("!");
  parameterExport->SetGuidance
  ("By default, will take a default value or the last \"/vis/ogl/set/printFilename\""
   " value if set.");
  fpCommandExport->SetParameter(parameterExport);
  parameterExport = new G4UIparameter ("width", 'd', omitable = true);
  parameterExport->SetGuidance
  ("By default, will take the current width of the viewer or \"/vis/ogl/set/printSize\""
   "\nif set. This parameter is only useful for eps/pdf/svg/ps formats !");
  parameterExport->SetDefaultValue(-1);
  fpCommandExport->SetParameter(parameterExport);
  parameterExport = new G4UIparameter ("height", 'd', omitable = true);
  parameterExport->SetGuidance
  ("By default, will take the current height of the viewer or \"/vis/ogl/set/printSize\""
   "\nif set. This parameter is only useful for eps/pdf/svg/ps formats !");
  parameterExport->SetDefaultValue(-1);
  fpCommandExport->SetParameter(parameterExport);

  fpCommandFlushAt = new G4UIcommand("/vis/ogl/flushAt", this);
  fpCommandFlushAt->SetGuidance
  ("Controls the rate at which graphics primitives are flushed to screen.");
  fpCommandFlushAt->SetGuidance
  ("Flushing to screen is an expensive operation so to speed drawing choose"
   "\nan action suitable for your application.  Note that detectors are flushed"
   "\nto screen anyway at end of drawing, and events are flushed to screen"
   "\nanyway depending on /vis/scene/endOfEventAction and endOfRunAction.");
  fpCommandFlushAt->SetGuidance
  ("For NthPrimitive and NthEvent the second parameter N is operative.");
  fpCommandFlushAt->SetGuidance
  ("For \"never\", detectors and events are still flushed as described above.");
  G4UIparameter* parameterFlushAt;
  parameterFlushAt = new G4UIparameter ("action", 's', omitable = true);
  parameterFlushAt->SetParameterCandidates
  ("endOfEvent endOfRun eachPrimitive NthPrimitive NthEvent never");
  parameterFlushAt->SetDefaultValue("NthEvent");
  fpCommandFlushAt->SetParameter(parameterFlushAt);
  parameterFlushAt = new G4UIparameter ("N", 'i', omitable = true);
  parameterFlushAt->SetDefaultValue(100);
  fpCommandFlushAt->SetParameter(parameterFlushAt);

  fpDirectorySet = new G4UIdirectory ("/vis/ogl/set/");
  fpDirectorySet->SetGuidance("G4OpenGLViewer set commands.");

  fpCommandDisplayListLimit =
    new G4UIcmdWithoutParameter("/vis/ogl/set/displayListLimit", this);
  fpCommandDisplayListLimit->SetGuidance
  ("This command is no longer relevant. There is no longer any limit on the"
   "\nnumber of display lists - except, of course, the available memory in"
   "\nyour computer. Keep an eye on that. Good luck!");

  fpCommandExportFormat =
  new G4UIcommand("/vis/ogl/set/exportFormat", this);
  fpCommandExportFormat->SetGuidance ("Set export format");
  fpCommandExportFormat->SetGuidance ("By default, pdf/eps/svg/ps are available. Depending of viewers several other format are available.");
  fpCommandExportFormat->SetGuidance ("Try /vis/ogl/set/exportFormat without parameters to see them.");
  fpCommandExportFormat->SetGuidance ("Changing format will reset the incremental suffix to 0.");
  G4UIparameter* parameterExportFormat;
  parameterExportFormat = new G4UIparameter ("format", 's', omitable = true);
  parameterExportFormat->SetDefaultValue("");
  fpCommandExportFormat->SetParameter(parameterExportFormat);

  fpCommandPrintFilename =
    new G4UIcommand("/vis/ogl/set/printFilename", this);
  fpCommandPrintFilename->SetGuidance ("Set print filename");
  fpCommandPrintFilename->SetGuidance ("Setting 'incremental' will increment filename by one at each new print, starting at 0");
  G4UIparameter* parameterPrintFilename;
  parameterPrintFilename = new G4UIparameter ("name", 's', omitable = true);
  parameterPrintFilename->SetDefaultValue("G4OpenGL");
  fpCommandPrintFilename->SetParameter(parameterPrintFilename);
  parameterPrintFilename = new G4UIparameter ("incremental", 'b', omitable = true);
  parameterPrintFilename->SetDefaultValue(1);
  fpCommandPrintFilename->SetParameter(parameterPrintFilename);
  
  fpCommandPrintMode = new G4UIcmdWithAString("/vis/ogl/set/printMode",this);
  fpCommandPrintMode->SetGuidance("Set print mode, only available for \"ps\" format");
  fpCommandPrintMode->SetParameterName("print_mode",omitable = true);
  fpCommandPrintMode->SetCandidates("vectored pixmap");
  fpCommandPrintMode->SetDefaultValue("vectored");

  fpCommandPrintSize =
    new G4UIcommand("/vis/ogl/set/printSize", this);
  fpCommandPrintSize->SetGuidance ("Set print size");
  fpCommandPrintSize->SetGuidance ("Tip : -1 will mean 'print size' = 'window size'");
  fpCommandPrintSize->SetGuidance ("       Setting size greater than your maximum graphic card capacity , will set the size to maximum  size.");
  G4UIparameter* parameterPrintSize;
  parameterPrintSize = new G4UIparameter ("width", 'd', omitable = false);
  parameterPrintSize->SetDefaultValue(-1);
  fpCommandPrintSize->SetParameter(parameterPrintSize);
  parameterPrintSize = new G4UIparameter ("height", 'd', omitable = false);
  parameterPrintSize->SetDefaultValue(-1);
  fpCommandPrintSize->SetParameter(parameterPrintSize);

  fpCommandTransparency =
    new G4UIcmdWithABool("/vis/ogl/set/transparency", this);
  fpCommandTransparency->SetGuidance
    ("True/false to enable/disable rendering of transparent objects.");
  fpCommandTransparency->SetParameterName
    ("transparency-enabled", omitable = true);
  fpCommandTransparency->SetDefaultValue(true);
}

G4OpenGLViewerMessenger::~G4OpenGLViewerMessenger ()
{
  delete fpCommandTransparency;
  delete fpCommandPrintSize;
  delete fpCommandPrintMode;
  delete fpCommandPrintFilename;
  delete fpCommandExportFormat;
  delete fpCommandDisplayListLimit;
  delete fpDirectorySet;
  delete fpCommandFlushAt;
  delete fpCommandExport;
  delete fpDirectory;

  delete fpInstance;
}

void G4OpenGLViewerMessenger::SetNewValue
(G4UIcommand* command, G4String newValue)
{
  G4VisManager* pVisManager = G4VisManager::GetInstance();

  G4VViewer* pViewer = pVisManager->GetCurrentViewer();
  if (!pViewer) {
    G4cout <<
      "G4OpenGLViewerMessenger::SetNewValue: No current viewer."
      "\n  \"/vis/open\", or similar, to get one."
           << G4endl;
    return;
  }

  G4VSceneHandler* pSceneHandler = pViewer->GetSceneHandler();
  if (!pSceneHandler) {
    G4cout <<
    "G4OpenGLViewerMessenger::SetNewValue: This viewer has no scene handler."
    "\n  Shouldn't happen - please report circumstances."
    "\n  (Viewer is \"" << pViewer->GetName() << "\".)"
    "\n  Try \"/vis/open\", or similar, to get one."
    << G4endl;
    return;
  }
  
  G4OpenGLViewer* pOGLViewer = dynamic_cast<G4OpenGLViewer*>(pViewer);
  if (!pOGLViewer) {
    G4cout <<
      "G4OpenGLViewerMessenger::SetNewValue: Current viewer is not of type"
      "\n  OGL.  (It is \""
	   << pViewer->GetName() <<
      "\".)\n  Use \"/vis/viewer/select\" or \"/vis/open\"."
           << G4endl;
    return;
  }

  G4OpenGLSceneHandler* pOGLSceneHandler =
  dynamic_cast<G4OpenGLSceneHandler*>(pSceneHandler);
  if (!pOGLSceneHandler) {
    G4cout <<
    "G4OpenGLViewerMessenger::SetNewValue: Current scene handler is not of type"
    "\n  OGL.  (Viewer is \"" << pViewer->GetName() << "\".)"
    "\n  (Scene handler is \"" << pSceneHandler->GetName() << "\".)"
    "\n  Use \"/vis/sceneHandler/list\" and \"/vis/sceneHandler/select\""
    "\n  or \"/vis/open\"."
    << G4endl;
    return;
  }

  if (command == fpCommandExport)
  {
    G4String name;
    G4int width,height;
    std::istringstream iss(newValue);
    iss >> name >> width >> height;
    pOGLViewer->exportImage(name, width, height);
    
    if (pOGLViewer->fVP.IsAutoRefresh())
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    return;
  }
  
  if (command == fpCommandExportFormat)
  {
    G4String name;
    std::istringstream iss(newValue);
    iss >> name;
    pOGLViewer->setExportImageFormat(name);

    return;
  }

  if (command == fpCommandFlushAt)
  {
    static G4bool firstTime = true;
    static std::map<G4String,G4OpenGLSceneHandler::FlushAction> actionMap;
    if (firstTime) {
      actionMap["endOfEvent"]    = G4OpenGLSceneHandler::endOfEvent;
      actionMap["endOfRun"]      = G4OpenGLSceneHandler::endOfRun;
      actionMap["eachPrimitive"] = G4OpenGLSceneHandler::eachPrimitive;
      actionMap["NthPrimitive"]  = G4OpenGLSceneHandler::NthPrimitive;
      actionMap["NthEvent"]      = G4OpenGLSceneHandler::NthEvent;
      actionMap["never"]         = G4OpenGLSceneHandler::never;
      firstTime = false;
    }
    G4String action;
    G4int entitiesFlushInterval;
    std::istringstream iss(newValue);
    iss >> action >> entitiesFlushInterval;
    pOGLSceneHandler->SetFlushAction(actionMap[action]);
    pOGLSceneHandler->SetEntitiesFlushInterval(entitiesFlushInterval);
    return;
  }

  if (command == fpCommandPrintFilename)
    {
      G4String name;
      G4bool inc;
      std::istringstream iss(newValue);
      iss >> name
    >> inc;
      pOGLViewer->setExportFilename(name,inc);
      return;
    }

  if (command == fpCommandPrintMode)
    {
      if (newValue == "vectored") pOGLViewer->fVectoredPs = true;

      if (newValue == "pixmap") pOGLViewer->fVectoredPs = false;
      return;
    }

  if (command == fpCommandPrintSize)
    {
      G4int width,height;
      std::istringstream iss(newValue);
      iss >> width
    >> height;
      pOGLViewer->setExportSize(width,height);
      return;
    }

  if (command == fpCommandTransparency)
    {
      pOGLViewer->transparency_enabled = command->ConvertToBool(newValue);
      if (pOGLViewer->fVP.IsAutoRefresh())
  G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
      return;
    }

  // Stored viewer commands
  G4OpenGLStoredViewer* pOGLSViewer =
    dynamic_cast<G4OpenGLStoredViewer*>(pViewer);

  if (!pOGLSViewer)
    {
      G4cout <<
  "G4OpenGLViewerMessenger::SetNewValue: Current viewer is not of type OGLS."
  "\n  (It is \"" << pViewer->GetName() << "\".)"
  "\n  This feature is only implemented for OGL Stored viewers."
  "\n  Use \"/vis/viewer/select\" or \"/vis/open OGLS...\"."
	     << G4endl;
      return;
    }

  // Scene handler commands
  G4OpenGLStoredSceneHandler* pOGLSSceneHandler =
    dynamic_cast<G4OpenGLStoredSceneHandler*>(pViewer->GetSceneHandler());

  if (!pOGLSSceneHandler) {
    G4cout <<
  "G4OpenGLViewerMessenger::SetNewValue: Current scene handler is not of type"
  "\n  OGLS (Stored).  (Viewer is \"" << pViewer->GetName() << "\".)"
  "\n  (Scene handler is \"" << pSceneHandler->GetName() << "\".)"
  "\n  This feature is only implemented for OGL Stored"
  "\n  scene handlers.  Use \"/vis/viewer/select\" or \"/vis/open OGLS...\"."
           << G4endl;
    return;
  }

  if (command == fpCommandDisplayListLimit)
    {
      G4cerr << command->GetGuidanceLine(0) << G4endl;
      return;
    }
}
