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
// $Id: NTSTDetectorMessenger.cc,v 1.2 2003-12-09 15:35:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "NTSTDetectorMessenger.hh"

#include "NTSTDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTDetectorMessenger::NTSTDetectorMessenger(NTSTDetectorConstruction * NTSTDet)
:NTSTDetector(NTSTDet)
{ 
  NTSTdetDir = new G4UIdirectory("/NTST/");
  NTSTdetDir->SetGuidance("NTST detector control.");

  InputFileNameCmd = new G4UIcmdWithAString("/NTST/setInputFile", this);
  InputFileNameCmd->SetGuidance("Set input file name");
  InputFileNameCmd->SetParameterName("File",true);
  InputFileNameCmd->SetDefaultValue("SVT.dat");

  DisableDet = new G4UIcmdWithAString("/NTST/disable", this);
  DisableDet->SetGuidance("disable detetector");
  DisableDet->SetCandidates("none SVT DCH all");
  DisableDet->SetDefaultValue("none");

  DebugCmd = new G4UIcmdWithAnInteger("/NTST/setDebug",this);
  DebugCmd->SetGuidance("Set debug flag.");
  DebugCmd->SetParameterName("Debug",true);

  NSubLayer = new G4UIcmdWithAnInteger("/NTST/setNSubLayer",this);
  NSubLayer->SetGuidance("Set the number of SVT sublayers.");
  NSubLayer->SetParameterName("NSubLay",true);
  NSubLayer->SetDefaultValue(7);
  NSubLayer->SetRange("NSubLay<8");

#if 0
  MinimumDriverStep
      = new G4UIcmdWithADoubleAndUnit("/NTST/setOuterRadius",this);
  MinimumDriverStep->SetGuidance("Set Minimum Step for ");
  MinimumDriverStep->SetParameterName("MinimumStep",false,false);
  MinimumDriverStep->SetDefaultValue(0.1);
  MinimumDriverStep->SetDefaultUnit("mm");
  MinimumDriverStep->SetRange("Radius>0.");
#endif
  
  MotherOuterRadius
      = new G4UIcmdWithADoubleAndUnit("/NTST/setOuterRadius",this);
  MotherOuterRadius->SetGuidance("Set outer radius of the SVT mother volume");
  MotherOuterRadius->SetParameterName("Radius",false,false);
  MotherOuterRadius->SetDefaultValue(19.);
  MotherOuterRadius->SetDefaultUnit("cm");
  MotherOuterRadius->SetRange("Radius>0.");
  
  fieldStat = new G4UIcmdWithoutParameter("/NTST/getFieldStats",this);
  fieldStat->SetGuidance( "Return number calls to field routine" );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTDetectorMessenger::~NTSTDetectorMessenger()
{
  delete DebugCmd; delete MotherOuterRadius;
  delete NTSTdetDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NTSTDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == DebugCmd )
   { NTSTDetector->SetDebugCmd(DebugCmd->GetNewIntValue(newValue));}
   
  if( command == MotherOuterRadius )
   { NTSTDetector->
         SetOuterRadius(MotherOuterRadius->GetNewDoubleValue(newValue));}

  if( command == NSubLayer )
   { NTSTDetector->
         SetNSubLayer(NSubLayer->GetNewIntValue(newValue));}

  if (command == fieldStat) NTSTDetector->GetFieldCallStats();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
