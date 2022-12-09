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
//
//
//


#include "G4RTMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4RTSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "G4VisManager.hh"
#include "G4RayTracerViewer.hh"
#include "G4TheRayTracer.hh"

#define G4warn G4cout

G4RTMessenger* G4RTMessenger::fpInstance = 0;

G4RTMessenger* G4RTMessenger::GetInstance
(G4TheRayTracer* p1)
{
  if (!fpInstance) fpInstance = new G4RTMessenger(p1);
  return fpInstance;
}

G4RTMessenger::G4RTMessenger(G4TheRayTracer* p1)
{
  theDefaultTracer = p1;
  theTracer = theDefaultTracer;

  rayDirectory = new G4UIdirectory("/vis/rayTracer/");
  rayDirectory->SetGuidance("RayTracer commands.");

  fileCmd = new G4UIcmdWithAString("/vis/rayTracer/trace",this);
  fileCmd->SetGuidance("Start the ray tracing.");
  fileCmd->SetGuidance("Define the name of output JPEG file.");
  fileCmd->SetParameterName("fileName",true);
  fileCmd->SetDefaultValue("g4RayTracer.jpeg");
  fileCmd->AvailableForStates(G4State_Idle);

  columnCmd = new G4UIcmdWithAnInteger("/vis/rayTracer/column",this);
  columnCmd->SetGuidance("Define the number of horizontal pixels.");
  columnCmd->SetParameterName("nPixel",false);
  columnCmd->SetRange("nPixel > 0");

  rowCmd = new G4UIcmdWithAnInteger("/vis/rayTracer/row",this);
  rowCmd->SetGuidance("Define the number of vertical pixels.");
  rowCmd->SetParameterName("nPixel",false);
  rowCmd->SetRange("nPixel > 0");

  targetCmd = new G4UIcmdWith3VectorAndUnit("/vis/rayTracer/target",this);
  targetCmd->SetGuidance("Define the center position of the target.");
  targetCmd->SetParameterName("X","Y","Z",true);
  targetCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  targetCmd->SetDefaultUnit("m");

  eyePosCmd = new G4UIcmdWith3VectorAndUnit("/vis/rayTracer/eyePosition",this);
  eyePosCmd->SetGuidance("Define the eye position.");
  eyePosCmd->SetGuidance("Eye direction is calculated from (target - eyePosition).");
  eyePosCmd->SetParameterName("X","Y","Z",true);
  eyePosCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  eyePosCmd->SetDefaultUnit("m");

  lightCmd = new G4UIcmdWith3Vector("/vis/rayTracer/lightDirection",this);
  lightCmd->SetGuidance("Define the direction of illumination light.");
  lightCmd->SetGuidance("The vector needs not to be a unit vector, but it must not be a zero vector.");
  lightCmd->SetParameterName("Px","Py","Pz",true);
  lightCmd->SetDefaultValue(G4ThreeVector(0.1,0.2,0.3));
  lightCmd->SetRange("Px != 0 || Py != 0 || Pz != 0");

  spanXCmd = new G4UIcmdWithADoubleAndUnit("/vis/rayTracer/span",this);
  spanXCmd->SetGuidance("Define the angle per 100 pixels.");
  spanXCmd->SetParameterName("span",true);
  spanXCmd->SetDefaultValue(50.);
  spanXCmd->SetDefaultUnit("deg");
  spanXCmd->SetRange("span>0.");

  headCmd = new G4UIcmdWithADoubleAndUnit("/vis/rayTracer/headAngle",this);
  headCmd->SetGuidance("Define the head direction.");
  headCmd->SetParameterName("headAngle",true);
  headCmd->SetDefaultValue(270.);
  headCmd->SetDefaultUnit("deg");
  headCmd->SetRange("headAngle>=0. && headAngle<360.");

  attCmd = new G4UIcmdWithADoubleAndUnit("/vis/rayTracer/attenuation",this);
  attCmd->SetGuidance("Define the attenuation length for transparent material.");
  attCmd->SetGuidance("Note that this value is independent to the attenuation length for the optical photon processes.");
  attCmd->SetParameterName("Length",true);
  attCmd->SetDefaultValue(1.0);
  attCmd->SetDefaultUnit("m");
  attCmd->SetRange("Length > 0.");

  distCmd = new G4UIcmdWithABool("/vis/rayTracer/distortion",this);
  distCmd->SetGuidance("Distortion effect of the fish eye lens.");
  distCmd->SetParameterName("flag",true);
  distCmd->SetDefaultValue(false);

  transCmd = new G4UIcmdWithABool("/vis/rayTracer/ignoreTransparency",this);
  transCmd->SetGuidance("Ignore transparency even if the alpha of G4Colour < 1.");
  transCmd->SetParameterName("flag",true);
  transCmd->SetDefaultValue(false);

  bkgColCmd = new G4UIcmdWith3Vector("/vis/rayTracer/backgroundColour",this);
  bkgColCmd->SetGuidance("Command has been deprecated.  Use /vis/viewer/set/background instead.");
  bkgColCmd->SetParameterName("red","green","blue",true);
  bkgColCmd->SetDefaultValue(G4ThreeVector(1.,1.,1.));
}

G4RTMessenger::~G4RTMessenger()
{
  delete columnCmd;
  delete rowCmd;
  delete targetCmd;
  delete eyePosCmd;
  delete lightCmd;
  delete spanXCmd;
  delete headCmd;
  delete attCmd;
  delete distCmd;
  delete transCmd;
  delete fileCmd;
  delete bkgColCmd;
  delete rayDirectory;
}

G4String G4RTMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String currentValue;
  if(command==columnCmd)
  { currentValue = columnCmd->ConvertToString(theTracer->GetNColumn()); }
  else if(command==rowCmd)
  { currentValue = rowCmd->ConvertToString(theTracer->GetNRow()); }
  else if(command==targetCmd)
  { currentValue = targetCmd->ConvertToString(theTracer->GetTargetPosition(),"m"); }
  else if(command==eyePosCmd)
  { currentValue = eyePosCmd->ConvertToString(theTracer->GetEyePosition(),"m"); }
  else if(command==lightCmd)
  { currentValue = lightCmd->ConvertToString(theTracer->GetLightDirection()); }
  else if(command==spanXCmd)
  { currentValue = spanXCmd->ConvertToString(theTracer->GetViewSpan(),"deg"); }
  else if(command==headCmd)
  { currentValue = headCmd->ConvertToString(theTracer->GetHeadAngle(),"deg"); }
  else if(command==attCmd)
  { currentValue = attCmd->ConvertToString(theTracer->GetAttenuationLength(),"m");}
  else if(command==distCmd)
  { currentValue = distCmd->ConvertToString(theTracer->GetDistortion()); }
  else if(command==transCmd)
  { currentValue = transCmd->ConvertToString(G4RTSteppingAction::GetIgnoreTransparency()); }
  else if(command==bkgColCmd)
  { currentValue = bkgColCmd->ConvertToString(theTracer->GetBackgroundColour()); }
  return currentValue;
}

void G4RTMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  G4VisManager* pVisManager = G4VisManager::GetInstance();

  theTracer = theDefaultTracer;

  G4VViewer* pVViewer = pVisManager->GetCurrentViewer();
  if (pVViewer) {
    G4RayTracerViewer* pViewer = dynamic_cast<G4RayTracerViewer*>(pVViewer);
    if (pViewer) {
      theTracer = pViewer->GetTracer();
    } else {
      G4warn <<
	"G4RTMessenger::SetNewValue: Current viewer is not of type RayTracer."
	"\n  Use \"/vis/viewer/select\" or \"/vis/open\"."
	     << G4endl;
    }
  }

  if (theTracer == theDefaultTracer) {
    G4warn <<
"G4RTMessenger::SetNewValue: No valid current viewer. Using default RayTracer."
	   << G4endl;
  }

  if(command==columnCmd)
  { theTracer->SetNColumn(columnCmd->GetNewIntValue(newValue)); }
  else if(command==rowCmd)
  { theTracer->SetNRow(rowCmd->GetNewIntValue(newValue)); }
  else if(command==targetCmd)
  { theTracer->SetTargetPosition(targetCmd->GetNew3VectorValue(newValue)); }
  else if(command==eyePosCmd)
  { theTracer->SetEyePosition(eyePosCmd->GetNew3VectorValue(newValue)); }
  else if(command==lightCmd)
  { theTracer->SetLightDirection(lightCmd->GetNew3VectorValue(newValue)); }
  else if(command==spanXCmd)
  { theTracer->SetViewSpan(spanXCmd->GetNewDoubleValue(newValue)); }
  else if(command==headCmd)
  { theTracer->SetHeadAngle(headCmd->GetNewDoubleValue(newValue)); }
  else if(command==attCmd)
  { theTracer->SetAttenuationLength(attCmd->GetNewDoubleValue(newValue)); }
  else if(command==distCmd)
  { theTracer->SetDistortion(distCmd->GetNewBoolValue(newValue)); }
  else if(command==bkgColCmd)
  {
	G4warn << "WARNING: /vis/rayTracer/backgroundColour has been deprecated."
	"\n  Use \"/vis/viewer/set/background\" instead."
		<< G4endl;
  }
  else if(command==transCmd)
  { G4RTSteppingAction::SetIgnoreTransparency(transCmd->GetNewBoolValue(newValue)); }
  else if(command==fileCmd)
  { theTracer->Trace(newValue); }
}


 


