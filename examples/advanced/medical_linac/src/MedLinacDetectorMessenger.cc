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
// $Id: MedLinacDetectorMessenger.cc,v 1.4 2005/07/03 23:27:37 mpiergen Exp $
//
//  Code developed by: M. Piergentili

#include "MedLinacDetectorMessenger.hh"
#include "MedLinacDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//*********************************************************************

MedLinacDetectorMessenger::MedLinacDetectorMessenger(
                                           MedLinacDetectorConstruction* MedLinacDet)
  :MedLinacDetector(MedLinacDet)
{ 

  G4cout <<"==================DetectorMessenger  "<<G4endl;
  MedLinacDir = new G4UIdirectory("/Jaws/");
  MedLinacDir->SetGuidance("jaws position");
  
  MedLinacDir = new G4UIdirectory("/Phantom/");
  MedLinacDir->SetGuidance("phantom parameters");


  X1Dir = new G4UIdirectory("/Jaws/X1/");
  X2Dir = new G4UIdirectory("/Jaws/X2/");
  Y1Dir = new G4UIdirectory("/Jaws/Y1/");
  Y2Dir = new G4UIdirectory("/Jaws/Y2/");


  JawX1PosCmd = new G4UIcmdWithADoubleAndUnit("/Jaws/X1/DistanceFromAxis",this);
  JawX1PosCmd->SetGuidance("Set jawsx1 position (negative)");
  JawX1PosCmd->SetParameterName("JawsX1Pos_x",false);
  JawX1PosCmd->SetRange("JawsX1Pos_x<=0. && JawsX1Pos_x>-1000.");
  JawX1PosCmd->SetDefaultUnit( "cm" );
  JawX1PosCmd->SetUnitCategory("Length");
  JawX1PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  JawX2PosCmd = new G4UIcmdWithADoubleAndUnit("/Jaws/X2/DistanceFromAxis",this);
  JawX2PosCmd->SetGuidance("Set jawsx2 position (positive)");
  JawX2PosCmd->SetParameterName("JawsX2Pos_x",false);
  JawX2PosCmd->SetRange("JawsX2Pos_x>=0. && JawsX2Pos_x<1000.");
  JawX2PosCmd->SetDefaultUnit( "cm" );
  JawX2PosCmd->SetUnitCategory("Length");
  JawX2PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  JawY1PosCmd = new G4UIcmdWithADoubleAndUnit("/Jaws/Y1/DistanceFromAxis",this);
  JawY1PosCmd->SetGuidance("Set jawsy1 position (negative)");
  JawY1PosCmd->SetParameterName("JawsY1Pos_y",false);
  JawY1PosCmd->SetRange("JawsY1Pos_y<=0. && JawsY1Pos_y>-1000.");
  JawY1PosCmd->SetDefaultUnit( "cm" );
  JawY1PosCmd->SetUnitCategory("Length");
  JawY1PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  JawY2PosCmd = new G4UIcmdWithADoubleAndUnit("/Jaws/Y2/DistanceFromAxis",this);
  JawY2PosCmd->SetGuidance("Set jawsy2 position (positive)");
  JawY2PosCmd->SetParameterName("JawsY2Pos_y",false);
  JawY2PosCmd->SetRange("JawsY2Pos_y>=0. && JawsY2Pos_y<1000.");
  JawY2PosCmd->SetDefaultUnit( "cm" );
  JawY2PosCmd->SetUnitCategory("Length");
  JawY2PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PhantomDimCmd = new G4UIcmdWithADoubleAndUnit("/Phantom/dimension",this);
  PhantomDimCmd->SetGuidance("Set phantom dimension (cm)");
  PhantomDimCmd->SetParameterName("phantomDim",false);
  PhantomDimCmd->SetDefaultUnit( "cm" );
  PhantomDimCmd->SetUnitCategory("Length");
  PhantomDimCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NVoxelsCmd = new G4UIcmdWithAnInteger("/Phantom/Nvoxels",this);
  NVoxelsCmd->SetGuidance("Set number of voxels along one axis");
  NVoxelsCmd->SetParameterName("numberOfVoxels",false);
  NVoxelsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MaxStepCmd = new G4UIcmdWithADoubleAndUnit("/Phantom/maxStep",this);
  MaxStepCmd->SetGuidance("Set max step in the phantom (mm)");
  MaxStepCmd->SetParameterName("maxStep",false);
  MaxStepCmd->SetDefaultUnit( "mm" );
  MaxStepCmd->SetUnitCategory("Length");
  MaxStepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/Jaws/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);


}

//****************************************************************************
MedLinacDetectorMessenger::~MedLinacDetectorMessenger()
{
  delete JawX1PosCmd;
  delete JawX2PosCmd;
  delete JawY1PosCmd;
  delete JawY2PosCmd;
  delete PhantomDimCmd;
  delete NVoxelsCmd;
  delete MaxStepCmd;
  delete UpdateCmd;

  delete MedLinacDir;
  
  delete X1Dir;
  delete X2Dir;
  delete Y1Dir;
  delete Y2Dir;
}

//****************************************************************************

void MedLinacDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == JawX1PosCmd )
   { MedLinacDetector->SetJawX1Pos_x(JawX1PosCmd->GetNewDoubleValue(newValue));}

  if( command == JawX2PosCmd )
   { MedLinacDetector->SetJawX2Pos_x(JawX2PosCmd->GetNewDoubleValue(newValue));}

  if( command == JawY1PosCmd )
   { MedLinacDetector->SetJawY1Pos_y(JawY1PosCmd->GetNewDoubleValue(newValue));}

  if( command == JawY2PosCmd )
   { MedLinacDetector->SetJawY2Pos_y(JawY2PosCmd->GetNewDoubleValue(newValue));}

  if( command == PhantomDimCmd )
   { MedLinacDetector->SetPhantomDim(PhantomDimCmd->GetNewDoubleValue(newValue));}

  if( command == NVoxelsCmd )
   { MedLinacDetector->SetNumberOfVoxels(NVoxelsCmd->GetNewIntValue(newValue));  }

  if( command == MaxStepCmd )
   { MedLinacDetector->SetMaxStep(MaxStepCmd->GetNewDoubleValue(newValue));}

  if( command == UpdateCmd )
   { MedLinacDetector->UpdateGeometry(); }

}
//****************************************************************************
