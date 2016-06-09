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
// $Id: G4OpenGLViewerMessenger.cc,v 1.8 2007/05/16 15:59:58 allison Exp $
// GEANT4 tag $Name: geant4-09-00 $

#include "G4OpenGLViewerMessenger.hh"

#include "G4OpenGLViewer.hh"
#include "G4OpenGLStoredViewer.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
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

  fpCommandPrintEPS =
    new G4UIcmdWithoutParameter("/vis/ogl/printEPS", this);
  fpCommandPrintEPS->SetGuidance("Print Encapsulated PostScript file.");
  fpCommandPrintEPS->SetGuidance
    ("Generates files with names G4OpenGL_n.eps, where n is a sequence"
     "\nnumber, starting at 0.");

  fpDirectorySet = new G4UIdirectory ("/vis/ogl/set/");
  fpDirectorySet->SetGuidance("G4OpenGLViewer set commands.");

  G4UIparameter* parameter;

  fpCommandDisplayHeadTime =
    new G4UIcommand("/vis/ogl/set/displayHeadTime", this);
  fpCommandDisplayHeadTime->SetGuidance
    ("Display head time of range in 2D text.");
  parameter = new G4UIparameter ("displayHeadTime", 'b', omitable = false);
  parameter->SetDefaultValue(false);
  fpCommandDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("screenX", 'd', omitable = true);
  parameter->SetGuidance("-1 < screenX < 1");
  parameter->SetParameterRange("screenX >= -1. && screenX <= 1.");
  parameter->SetDefaultValue(-0.9);
  fpCommandDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("screenY", 'd', omitable = true);
  parameter->SetGuidance("-1 < screenY < 1");
  parameter->SetParameterRange("screenY >= -1. && screenY <= 1.");
  parameter->SetDefaultValue(-0.9);
  fpCommandDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("screenSize", 'd', omitable = true);
  parameter->SetDefaultValue(24.);
  fpCommandDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("red", 'd', omitable = true);
  parameter->SetParameterRange("red >= 0. && red <= 1.");
  parameter->SetDefaultValue(0.);
  fpCommandDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetParameterRange("green >= 0. && green <= 1.");
  parameter->SetDefaultValue(1.);
  fpCommandDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetParameterRange("blue >= 0. && blue <= 1.");
  parameter->SetDefaultValue(1.);
  fpCommandDisplayHeadTime->SetParameter(parameter);

  fpCommandDisplayLightFront =
    new G4UIcommand("/vis/ogl/set/displayLightFront", this);
  fpCommandDisplayLightFront->SetGuidance
    ("Display the light front at head time.");
  fpCommandDisplayLightFront->SetGuidance
    ("Tip: The trajectories can appear of jump ahead of the light front"
     "\nbecause their time range overlaps the viewer's time range.  To"
     "\naverage out this discrete time effect, advance the light front by"
     "\nhalf the trajectories interval. E.g., if the trajectory time slice"
     "\ninterval is 0.01 ns:"
     "\n  /vis/ogl/set/displayLightFront true -90 0 0 mm -0.005 ns"
     "\nTo prevent them beating the light front at all:"
     "\n  /vis/ogl/set/displayLightFront true -90 0 0 mm -0.01 ns");
  parameter = new G4UIparameter ("displayLightFront", 'b', omitable = false);
  parameter->SetDefaultValue(false);
  fpCommandDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("originX", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  fpCommandDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("originY", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  fpCommandDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("originZ", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  fpCommandDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("space_unit", 's', omitable = true);
  parameter->SetDefaultValue("m");
  fpCommandDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("originT", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  fpCommandDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("time_unit", 's', omitable = true);
  parameter->SetDefaultValue("s");
  fpCommandDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("red", 'd', omitable = true);
  parameter->SetParameterRange("red >= 0. && red <= 1.");
  parameter->SetDefaultValue(0.);
  fpCommandDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetParameterRange("green >= 0. && green <= 1.");
  parameter->SetDefaultValue(1.);
  fpCommandDisplayLightFront->SetParameter(parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetParameterRange("blue >= 0. && blue <= 1.");
  parameter->SetDefaultValue(0.);
  fpCommandDisplayLightFront->SetParameter(parameter);

  fpCommandEndTime =
    new G4UIcommand("/vis/ogl/set/endTime", this);
  fpCommandEndTime->SetGuidance("Set end and range of track time.");
  parameter = new G4UIparameter ("end-time", 'd', omitable = false);
  parameter->SetDefaultValue(G4OPENGL_DBL_MAX);
  fpCommandEndTime->SetParameter(parameter);
  parameter = new G4UIparameter ("end-time-unit", 's', omitable = false);
  parameter->SetDefaultValue("ns");
  fpCommandEndTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range", 'd', omitable = true);
  parameter->SetDefaultValue(-1.);
  fpCommandEndTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range-unit", 's', omitable = true);
  parameter->SetDefaultValue("ns");
  fpCommandEndTime->SetParameter(parameter);

  fpCommandFade = new G4UIcmdWithADouble("/vis/ogl/set/fade", this);
  fpCommandFade->SetGuidance
    ("0: no fade; 1: maximum fade with time within range.");
  fpCommandFade->SetParameterName("fadefactor", omitable = false);
  fpCommandFade->SetRange("fadefactor>=0.&&fadefactor<=1.");
  fpCommandFade->SetDefaultValue(0.);

  fpCommandPrintMode = new G4UIcmdWithAString
    ("/vis/ogl/set/printMode",this);
  fpCommandPrintMode->SetGuidance("Set print mode");
  fpCommandPrintMode->SetParameterName("print_mode",omitable = true);
  fpCommandPrintMode->SetCandidates("vectored pixmap");
  fpCommandPrintMode->SetDefaultValue("vectored");

  fpCommandStartTime =
    new G4UIcommand("/vis/ogl/set/startTime", this);
  fpCommandStartTime->SetGuidance("Set start and range of track time.");
  parameter = new G4UIparameter ("start-time", 'd', omitable = false);
  parameter->SetDefaultValue(-G4OPENGL_DBL_MAX);
  fpCommandStartTime->SetParameter(parameter);
  parameter = new G4UIparameter ("start-time-unit", 's', omitable = false);
  parameter->SetDefaultValue("ns");
  fpCommandStartTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range", 'd', omitable = true);
  parameter->SetDefaultValue(-1.);
  fpCommandStartTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range-unit", 's', omitable = true);
  parameter->SetDefaultValue("ns");
  fpCommandStartTime->SetParameter(parameter);

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
  delete fpCommandPrintMode;
  delete fpCommandTransparency;
  delete fpCommandStartTime;
  delete fpCommandFade;
  delete fpCommandEndTime;
  delete fpCommandDisplayLightFront;
  delete fpCommandDisplayHeadTime;
  delete fpDirectorySet;
  delete fpCommandPrintEPS;
  delete fpDirectory;
}

void G4OpenGLViewerMessenger::SetNewValue
(G4UIcommand* command, G4String newValue)
{
  G4VisManager* pVisManager = G4VisManager::GetInstance();

  G4VViewer* pVViewer = pVisManager->GetCurrentViewer();

  if (!pVViewer) {
    G4cout <<
      "G4OpenGLViewerMessenger::SetNewValue: No current viewer."
      "\n  \"/vis/open\", or similar, to get one."
           << G4endl;
    return;
  }

  G4OpenGLViewer* pOGLViewer = dynamic_cast<G4OpenGLViewer*>(pVViewer);

  if (!pOGLViewer) {
    G4cout <<
      "G4OpenGLViewerMessenger::SetNewValue: Current viewer is not of type"
      "\n  OGL.  Use \"/vis/viewer/select\" or \"/vis/open\"."
           << G4endl;
    return;
  }

  if (command == fpCommandPrintEPS) 
    {
      // Keep copy of print_string to preserve Xm behaviour...
      char* tmp_string = new char[50];
      strcpy (tmp_string, pOGLViewer->print_string);
      // Make new print string...
      static G4int file_count = 0;
      std::ostringstream oss;
      oss << "G4OpenGL_" << file_count++ << ".eps";
      strcpy (pOGLViewer->print_string, oss.str().c_str());
      // Print eps file...
      pOGLViewer->print();
      // Restore print_string for Xm...
      strcpy (pOGLViewer->print_string, tmp_string);
      delete tmp_string;
      return;
    }

  G4OpenGLStoredViewer* pViewer =
    dynamic_cast<G4OpenGLStoredViewer*>(pVViewer);

  if (!pViewer) {
    G4cout <<
  "G4OpenGLViewerMessenger::SetNewValue: Current viewer is not of type OGLS."
  "\n  The time slice viewing feature is only implemented for OGL Stored"
  "\n  viewers at present.  Use \"/vis/viewer/select\" or \"/vis/open\"."
           << G4endl;
    return;
  }

  if (command == fpCommandDisplayHeadTime)
    {
      G4String display;
      G4double screenX, screenY, screenSize, red, green, blue;
      std::istringstream iss(newValue);
      iss >> display >> screenX >> screenY
	  >> screenSize >> red >> green >> blue;
      pViewer->fDisplayHeadTime = command->ConvertToBool(display);
      pViewer->fDisplayHeadTimeX = screenX;
      pViewer->fDisplayHeadTimeY = screenY;
      pViewer->fDisplayHeadTimeSize = screenSize;
      pViewer->fDisplayHeadTimeRed = red;
      pViewer->fDisplayHeadTimeGreen = green;
      pViewer->fDisplayHeadTimeBlue = blue;
    }

  if (command == fpCommandDisplayLightFront)
    {
      G4String display, originX, originY, originZ, unitS, originT, unitT;
      G4double red, green, blue;
      std::istringstream iss(newValue);
      iss >> display
	  >> originX >> originY >> originZ >> unitS
	  >> originT >> unitT
	  >> red >> green >> blue;
      pViewer->fDisplayLightFront = command->ConvertToBool(display);
      pViewer->fDisplayLightFrontX =
	command->ConvertToDimensionedDouble(G4String(originX + ' ' + unitS));
      pViewer->fDisplayLightFrontY =
	command->ConvertToDimensionedDouble(G4String(originY + ' ' + unitS));
      pViewer->fDisplayLightFrontZ =
	command->ConvertToDimensionedDouble(G4String(originZ + ' ' + unitS));
      pViewer->fDisplayLightFrontT =
	command->ConvertToDimensionedDouble(G4String(originT + ' ' + unitT));
      pViewer->fDisplayLightFrontRed = red;
      pViewer->fDisplayLightFrontGreen = green;
      pViewer->fDisplayLightFrontBlue = blue;
    }

  if (command == fpCommandEndTime)
    {
      G4String end_time_string, end_time_unit,
	time_range_string, time_range_unit;
      std::istringstream iss(newValue);
      iss >> end_time_string >> end_time_unit
	  >> time_range_string >> time_range_unit;
      pViewer->fEndTime = command->ConvertToDimensionedDouble
	(G4String(end_time_string + ' ' + end_time_unit));
      G4double timeRange = command->ConvertToDimensionedDouble
	(G4String(time_range_string + ' ' + time_range_unit));
      if (timeRange > 0.) {
	pViewer->fStartTime = pViewer->fEndTime - timeRange;
      }
      if (pViewer->fVP.IsAutoRefresh())
	G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    }

  if (command == fpCommandFade)
    {
      pViewer->fFadeFactor = command->ConvertToDouble(newValue);
      if (pViewer->fVP.IsAutoRefresh())
	G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    }

  if (command == fpCommandPrintMode)
    {
      if (newValue == "vectored") pViewer->vectored_ps = true;
      if (newValue == "pixmap") {
	pViewer->vectored_ps = false;
	if (pVisManager->GetVerbosity() >= G4VisManager::warnings) {
	  G4cout <<
	    "WARNING: Only implemented for X Windows at present."
		 << G4endl;
	}
      }
    }

  if (command == fpCommandStartTime)
    {
      G4String start_time_string, start_time_unit,
	time_range_string, time_range_unit;
      std::istringstream iss(newValue);
      iss >> start_time_string >> start_time_unit
	  >> time_range_string >> time_range_unit;
      pViewer->fStartTime = command->ConvertToDimensionedDouble
	(G4String(start_time_string + ' ' + start_time_unit));
      G4double timeRange = command->ConvertToDimensionedDouble
	(G4String(time_range_string + ' ' + time_range_unit));
      if (timeRange > 0.) {
	pViewer->fEndTime = pViewer->fStartTime + timeRange;
      }
      if (pViewer->fVP.IsAutoRefresh())
	G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    }

  if (command == fpCommandTransparency)
    {
      pViewer->transparency_enabled = command->ConvertToBool(newValue);
      if (pViewer->fVP.IsAutoRefresh())
	G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
    }

}
