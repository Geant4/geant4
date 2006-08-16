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
// $Id: G4OpenGLViewerMessenger.cc,v 1.2 2006-08-16 10:33:54 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4OpenGLViewerMessenger.hh"

#include "G4OpenGLViewer.hh"
#include "G4OpenGLStoredViewer.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
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
  parameter = new G4UIparameter ("screen-x", 'd', omitable = true);
  parameter->SetGuidance("-1 < screen-x < 1");
  parameter->SetDefaultValue(-0.9);
  fpCommandDisplayHeadTime->SetParameter(parameter);
  parameter = new G4UIparameter ("screen-y", 'd', omitable = true);
  parameter->SetGuidance("-1 < screen-y < 1");
  parameter->SetDefaultValue(-0.9);
  fpCommandDisplayHeadTime->SetParameter(parameter);

  fpCommandEndTime =
    new G4UIcmdWithADoubleAndUnit("/vis/ogl/set/endTime", this);
  fpCommandEndTime->SetGuidance("Set end of range of track time.");
  fpCommandEndTime->SetParameterName("end-time", omitable = false);
  fpCommandEndTime->SetDefaultValue(DBL_MAX);

  fpCommandFade = new G4UIcmdWithADouble("/vis/ogl/set/fade", this);
  fpCommandFade->SetGuidance
    ("0: no fade; 1: maximum fade with time within range.");
  fpCommandFade->SetParameterName("fadefactor", omitable = false);
  fpCommandFade->SetRange("fadefactor>=0.&&fadefactor<=1.");
  fpCommandFade->SetDefaultValue(0.);

  fpCommandStartTime =
    new G4UIcommand("/vis/ogl/set/startTime", this);
  fpCommandStartTime->SetGuidance("Set start and range of track time.");
  parameter = new G4UIparameter ("start-time", 'd', omitable = false);
  parameter->SetDefaultValue(-DBL_MAX);
  fpCommandStartTime->SetParameter(parameter);
  parameter = new G4UIparameter ("start-time-unit", 's', omitable = false);
  parameter->SetDefaultValue("ns");
  fpCommandStartTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range", 'd', omitable = true);
  parameter->SetDefaultValue(DBL_MAX);
  fpCommandStartTime->SetParameter(parameter);
  parameter = new G4UIparameter ("time-range-unit", 's', omitable = true);
  parameter->SetDefaultValue("ns");
  fpCommandStartTime->SetParameter(parameter);
}

G4OpenGLViewerMessenger::~G4OpenGLViewerMessenger ()
{
  delete fpCommandStartTime;
  delete fpCommandEndTime;
  delete fpDirectorySet;
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
      G4double screenX, screenY;
      std::istringstream iss(newValue);
      iss >> display >> screenX >> screenY;
      pViewer->fDisplayHeadTime = command->ConvertToBool(display);
      pViewer->fDisplayHeadTimeX = screenX;
      pViewer->fDisplayHeadTimeY = screenY;
    }

  if (command == fpCommandEndTime)
    {
      pViewer->fEndTime = command->ConvertToDimensionedDouble(newValue);
    }

  if (command == fpCommandFade)
    {
      pViewer->fFadeFactor = command->ConvertToDouble(newValue);
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
      pViewer->fEndTime = pViewer->fStartTime + timeRange;
   }

  G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");

}
