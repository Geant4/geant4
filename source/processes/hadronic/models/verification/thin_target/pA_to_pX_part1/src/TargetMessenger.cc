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
// $Id: TargetMessenger.cc,v 1.1 2003-05-27 13:44:49 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "TargetMessenger.hh"

#include "TargetConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"


TargetMessenger::TargetMessenger(TargetConstruction* Tgt)
:Target(Tgt)
{ 
  TargDir = new G4UIdirectory("/target/");
  TargDir->SetGuidance("Target detector control.");
      
  TgtMaterCmd = new G4UIcmdWithAString("/target/setTgtMat",this);
  TgtMaterCmd->SetGuidance("Select Material of the Absorber.");
  TgtMaterCmd->SetParameterName("choice",false);
  TgtMaterCmd->AvailableForStates(G4State_Idle);
  
  TgtThickCmd = new G4UIcmdWithADoubleAndUnit("/target/setTgtThick",this);
  TgtThickCmd->SetGuidance("Set Target Thickness");
  TgtThickCmd->SetParameterName("Size",false);
  TgtThickCmd->SetRange("Size>=0.");
  TgtThickCmd->SetUnitCategory("Length");
  TgtThickCmd->AvailableForStates(G4State_Idle);
  
  TgtRadiusCmd = new G4UIcmdWithADoubleAndUnit("/target/setTgtRadius",this);
  TgtRadiusCmd->SetGuidance("Set Target Radius");
  TgtRadiusCmd->SetParameterName("Size",false);
  TgtRadiusCmd->SetRange("Size>=0.");
  TgtRadiusCmd->SetUnitCategory("Length");  
  TgtRadiusCmd->AvailableForStates(G4State_Idle);
  
  NumEvtCmd = new G4UIcmdWithAnInteger("/target/setNumEvts",this);
  NumEvtCmd->SetGuidance("Set Number of Events to Be Run");
  //  NumEvtCmd->SetParameterName("Size",false);
  //  NumEvtCmd->SetRange("Size>=0.");
  NumEvtCmd->AvailableForStates(G4State_Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/target/update",this);
  UpdateCmd->SetGuidance("Update target geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/target/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_Idle);  
}


TargetMessenger::~TargetMessenger()
{
  delete TgtMaterCmd; 
  delete TgtThickCmd; 
  delete TgtRadiusCmd;
  delete NumEvtCmd;
  delete MagFieldCmd;
  delete TargDir;
}


void TargetMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == TgtMaterCmd ) { 
    Target->SetTargetMaterial(newValue);
  }
   
  if( command == TgtThickCmd ) { 
    Target->SetTargetThickness(TgtThickCmd->GetNewDoubleValue(newValue));
  }
   
  if( command == NumEvtCmd ) { 
    Target->SetNumRequestedEvents(NumEvtCmd->GetNewIntValue(newValue));
  }
  
  if( command == TgtRadiusCmd ) { 
    Target->SetTargetRadius(TgtRadiusCmd->GetNewDoubleValue(newValue));
  }
  
  if( command == UpdateCmd ) { 
    Target->UpdateGeometry(); 
  }

  if( command == MagFieldCmd ) { 
    Target->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));
  }
}








