// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01FieldMessenger.cc,v 1.3 2001-03-29 16:29:50 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "F01FieldMessenger.hh"

#include "F01ElectroMagneticField.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//////////////////////////////////////////////////////////////////////////////

F01FieldMessenger::F01FieldMessenger(F01ElectroMagneticField* pEMfield)
  :fEMfield(pEMfield)
{ 
  F01detDir = new G4UIdirectory("/field/");
  F01detDir->SetGuidance("F01 field tracking control.");

  StepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType",this);
  StepperCmd->SetGuidance("Select stepper type for magnetic field");
  StepperCmd->SetParameterName("choice",true);
  StepperCmd->SetDefaultValue(4);
  StepperCmd->AvailableForStates(PreInit,Idle);

 
  UpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/field/setFieldZ",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false,false);
  MagFieldCmd->SetDefaultUnit("tesla");
  MagFieldCmd->AvailableForStates(Idle); 
 
  MinStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep",this);  
  MinStepCmd->SetGuidance("Define minimal step");
  MinStepCmd->SetGuidance("Magnetic field will be in Z direction.");
  MinStepCmd->SetParameterName("min step",false,false);
  MinStepCmd->SetDefaultUnit("mm");
  MinStepCmd->AvailableForStates(Idle);  
       
  AbsMaterCmd = new G4UIcmdWithAString("/field/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",true);
  AbsMaterCmd->SetDefaultValue("Xe");
  AbsMaterCmd->AvailableForStates(Idle);


}

///////////////////////////////////////////////////////////////////////////////

F01FieldMessenger::~F01FieldMessenger()
{
  delete StepperCmd;
  delete MagFieldCmd;
  delete MinStepCmd;
  delete F01detDir;
  delete UpdateCmd;

  delete AbsMaterCmd; 
}

////////////////////////////////////////////////////////////////////////////
//
//

void F01FieldMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{ 
  if( command == StepperCmd )
  { 
    fEMfield->SetStepperType(StepperCmd->GetNewIntValue(newValue));
  }  
  if( command == UpdateCmd )
  { 
    fEMfield->UpdateField(); 
  }
  if( command == MagFieldCmd )
  { 
    fEMfield->SetFieldValue(MagFieldCmd->GetNewDoubleValue(newValue));
  }
  if( command == MinStepCmd )
  { 
    fEMfield->SetMinStep(MinStepCmd->GetNewDoubleValue(newValue));
  }
}

//
//
/////////////////////////////////////////////////////////////////////////



























