
#include "G4RTMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4RayTracer.hh"
#include "G4RTSteppingAction.hh"
#include "G4ThreeVector.hh"

G4RTMessenger::G4RTMessenger(G4RayTracer* p1,G4RTSteppingAction* p2)
{
  theTracer = p1;
  theSteppingAction = p2;

  rayDirectory = new G4UIdirectory("/vis/rayTracer/");
  rayDirectory->SetGuidance("RayTracer commands.");

  fileCmd = new G4UIcmdWithAString("/vis/rayTracer/trace",this);
  fileCmd->SetGuidance("Start the ray tracing.");
  fileCmd->SetGuidance("Define the name of output JPEG file.");
  fileCmd->SetParameterName("fileName",true);
  fileCmd->SetDefaultValue("g4RayTracer.jpeg");

  columnCmd = new G4UIcmdWithAnInteger("/vis/rayTracer/column",this);
  columnCmd->SetGuidance("Define the number of horizontal pixels.");
  columnCmd->SetParameterName("nPixel",false);
  columnCmd->SetRange("nPixel > 0");

  rowCmd = new G4UIcmdWithAnInteger("/vis/rayTracer/row",this);
  rowCmd->SetGuidance("Define the number of virtical pixels.");
  rowCmd->SetParameterName("nPixel",false);
  rowCmd->SetRange("nPixel > 0");

  targetCmd = new G4UIcmdWith3VectorAndUnit("/vis/rayTracer/target",this);
  targetCmd->SetGuidance("Define the center position of the target.");
  targetCmd->SetParameterName("X","Y","Z",true);
  targetCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  targetCmd->SetDefaultUnit("m");

  eyePosCmd = new G4UIcmdWith3VectorAndUnit("/vis/rayTracer/eyePosition",this);
  eyePosCmd->SetGuidance("Define the eye position.");
  eyePosCmd->SetGuidance("Eye direction is calsurated from (target - eyePosition).");
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
  headCmd->SetDefaultValue(0.);
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
  { currentValue = transCmd->ConvertToString(theSteppingAction->GetIgnoreTransparency()); }
  return currentValue;
}

void G4RTMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
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
  else if(command==transCmd)
  { theSteppingAction->SetIgnoreTransparency(transCmd->GetNewBoolValue(newValue)); }
  else if(command==fileCmd)
  { theTracer->Trace(newValue); }
}


 


