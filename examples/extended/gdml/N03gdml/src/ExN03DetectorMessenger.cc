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
// $Id: ExN03DetectorMessenger.cc,v 1.2 2002-12-05 01:06:59 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03DetectorMessenger.hh"

#include "ExN03DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03DetectorMessenger::ExN03DetectorMessenger(
                                           ExN03DetectorConstruction* ExN03Det)
:ExN03Detector(ExN03Det)
{ 
  GDMLDir = new G4UIdirectory("/gdml/");
  GDMLDir->SetGuidance("GDML control.");
  
  GDMLFile = new G4UIcmdWithAString("/gdml/input",this);
  GDMLFile->SetGuidance("Select GDML input");
  GDMLFile->SetParameterName("file",false);
  GDMLFile->AvailableForStates(G4State_Idle);

  GDMLSetup = new G4UIcmdWithAString("/gdml/setup",this);
  GDMLSetup->SetGuidance("Select geometry setup in GDML input");
  GDMLSetup->SetGuidance("If none is specified the GDML processor will choose the first found");
  GDMLSetup->SetParameterName("id",true);
  GDMLSetup->SetDefaultValue( "default" );
  GDMLSetup->AvailableForStates(G4State_Idle);

  GDMLVersion = new G4UIcmdWithAString("/gdml/version",this);
  GDMLVersion->SetGuidance("Select a particular version of a given geometry setup in GDML input");
  GDMLVersion->SetGuidance("If none is specified the GDML processor will choose the first found");
  GDMLVersion->SetParameterName("id",false);
  GDMLVersion->SetDefaultValue( "default" );
  GDMLVersion->AvailableForStates(G4State_Idle);
  
  ExN03detDir = new G4UIdirectory("/calor/");
  ExN03detDir->SetGuidance("ExN03 detector control.");
      
  AbsMaterCmd = new G4UIcmdWithAString("/calor/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",false);
  AbsMaterCmd->AvailableForStates(G4State_Idle);
  
  GapMaterCmd = new G4UIcmdWithAString("/calor/setGapMat",this);
  GapMaterCmd->SetGuidance("Select Material of the Gap.");
  GapMaterCmd->SetParameterName("choice",false);
  GapMaterCmd->AvailableForStates(G4State_Idle);
    
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("Size",false);
  AbsThickCmd->SetRange("Size>=0.");
  AbsThickCmd->SetUnitCategory("Length");
  AbsThickCmd->AvailableForStates(G4State_Idle);
  
  GapThickCmd = new G4UIcmdWithADoubleAndUnit("/calor/setGapThick",this);
  GapThickCmd->SetGuidance("Set Thickness of the Gap");
  GapThickCmd->SetParameterName("Size",false);
  GapThickCmd->SetRange("Size>=0.");
  GapThickCmd->SetUnitCategory("Length");  
  GapThickCmd->AvailableForStates(G4State_Idle);
  
  SizeYZCmd = new G4UIcmdWithADoubleAndUnit("/calor/setSizeYZ",this);
  SizeYZCmd->SetGuidance("Set tranverse size of the calorimeter");
  SizeYZCmd->SetParameterName("Size",false);
  SizeYZCmd->SetRange("Size>0.");
  SizeYZCmd->SetUnitCategory("Length");    
  SizeYZCmd->AvailableForStates(G4State_Idle);
  
  NbLayersCmd = new G4UIcmdWithAnInteger("/calor/setNbOfLayers",this);
  NbLayersCmd->SetGuidance("Set number of layers.");
  NbLayersCmd->SetParameterName("NbLayers",false);
  NbLayersCmd->SetRange("NbLayers>0 && NbLayers<500");
  NbLayersCmd->AvailableForStates(G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/calor/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/calor/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03DetectorMessenger::~ExN03DetectorMessenger()
{
  delete NbLayersCmd;
  delete AbsMaterCmd; delete GapMaterCmd;
  delete AbsThickCmd; delete GapThickCmd;
  delete SizeYZCmd;   delete UpdateCmd;
  delete MagFieldCmd;
  delete ExN03detDir;
  delete GDMLFile;
  delete GDMLSetup;
  delete GDMLVersion;
  delete GDMLDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  /*
  if( command == AbsMaterCmd )
   { ExN03Detector->SetAbsorberMaterial(newValue);}
   
  if( command == GapMaterCmd )
   { ExN03Detector->SetGapMaterial(newValue);}
  
  if( command == AbsThickCmd )
   { ExN03Detector->SetAbsorberThickness(AbsThickCmd
                                               ->GetNewDoubleValue(newValue));}
   
  if( command == GapThickCmd )
   { ExN03Detector->SetGapThickness(GapThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == SizeYZCmd )
   { ExN03Detector->SetCalorSizeYZ(SizeYZCmd->GetNewDoubleValue(newValue));}
   
  if( command == NbLayersCmd )
   { ExN03Detector->SetNbOfLayers(NbLayersCmd->GetNewIntValue(newValue));}
  
  if( command == UpdateCmd )
   { ExN03Detector->UpdateGeometry(); }

  if( command == MagFieldCmd )
   { ExN03Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
