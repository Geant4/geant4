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
// $Id: DetectorMessenger.cc,v 1.3 2006-06-29 16:58:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{ 
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance(" detector control.");
  
  detDir = new G4UIdirectory("/testem/det/");
  detDir->SetGuidance("detector construction commands");
      
  MaterCmd = new G4UIcmdWithAString("/testem/det/setMat",this);
  MaterCmd->SetGuidance("Select material of the box.");
  MaterCmd->SetParameterName("choice",false);
  MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SizeXCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeX",this);
  SizeXCmd->SetGuidance("Set sizeX of the absorber");
  SizeXCmd->SetParameterName("SizeX",false);
  SizeXCmd->SetRange("SizeX>0.");
  SizeXCmd->SetUnitCategory("Length");
  SizeXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SizeYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeYZ",this);
  SizeYZCmd->SetGuidance("Set sizeYZ of the absorber");
  SizeYZCmd->SetParameterName("SizeYZ",false);
  SizeYZCmd->SetRange("SizeYZ>0.");
  SizeYZCmd->SetUnitCategory("Length");
  SizeYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
        
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  TalNbCmd = new G4UIcmdWithAnInteger("/testem/det/tallyNumber",this);
  TalNbCmd->SetGuidance("Set number of Tallies.");
  TalNbCmd->SetParameterName("tallyNb",false);
  TalNbCmd->SetRange("tallyNb>=0");
  TalNbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TalDefCmd = new G4UIcommand("/testem/det/tallyDefinition",this);
  TalDefCmd->SetGuidance("Set tally nb, material, box dimensions.");
  TalDefCmd->SetGuidance("  tally number : from 1 to tallyNumber");
  TalDefCmd->SetGuidance("  material name");
  TalDefCmd->SetGuidance("  dimensions (3-vector with unit)");
  //
  G4UIparameter* TalNbPrm = new G4UIparameter("tallyNb",'i',false);
  TalNbPrm->SetGuidance("tally number : from 1 to tallyNumber");
  TalNbPrm->SetParameterRange("tallyNb>0");
  TalDefCmd->SetParameter(TalNbPrm);
  //
  G4UIparameter* MatPrm = new G4UIparameter("material",'s',false);
  MatPrm->SetGuidance("material name");
  TalDefCmd->SetParameter(MatPrm);
  //    
  G4UIparameter* SizeXPrm = new G4UIparameter("sizeX",'d',false);
  SizeXPrm->SetGuidance("sizeX");
  SizeXPrm->SetParameterRange("sizeX>0.");
  TalDefCmd->SetParameter(SizeXPrm);
  //    
  G4UIparameter* SizeYPrm = new G4UIparameter("sizeY",'d',false);
  SizeYPrm->SetGuidance("sizeY");
  SizeYPrm->SetParameterRange("sizeY>0.");
  TalDefCmd->SetParameter(SizeYPrm);
  //    
  G4UIparameter* SizeZPrm = new G4UIparameter("sizeZ",'d',false);
  SizeZPrm->SetGuidance("sizeZ");
  SizeZPrm->SetParameterRange("sizeZ>0.");
  TalDefCmd->SetParameter(SizeZPrm);    
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of dimensions");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitList);
  TalDefCmd->SetParameter(unitPrm);
  //
  TalDefCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TalPosiCmd = new G4UIcommand("/testem/det/tallyPosition",this);
  TalPosiCmd->SetGuidance("Set tally nb, position");
  TalPosiCmd->SetGuidance("  tally number : from 1 to tallyNumber");
  TalPosiCmd->SetGuidance("  position (3-vector with unit)");
  //
  G4UIparameter* TalNumPrm = new G4UIparameter("tallyNum",'i',false);
  TalNumPrm->SetGuidance("tally number : from 1 to tallyNumber");
  TalNumPrm->SetParameterRange("tallyNum>0");
  TalPosiCmd->SetParameter(TalNumPrm);
  //    
  G4UIparameter* PosiXPrm = new G4UIparameter("posiX",'d',false);
  PosiXPrm->SetGuidance("position X");
  TalPosiCmd->SetParameter(PosiXPrm);
  //    
  G4UIparameter* PosiYPrm = new G4UIparameter("posiY",'d',false);
  PosiYPrm->SetGuidance("position Y");
  TalPosiCmd->SetParameter(PosiYPrm);
  //
  G4UIparameter* PosiZPrm = new G4UIparameter("posiZ",'d',false);
  PosiZPrm->SetGuidance("position Z");
  TalPosiCmd->SetParameter(PosiZPrm);      
  //
  G4UIparameter* unitPr = new G4UIparameter("unit",'s',false);
  unitPr->SetGuidance("unit of position");
  unitPr->SetParameterCandidates(unitList);
  TalPosiCmd->SetParameter(unitPr);
  //
  TalPosiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
        
  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete MaterCmd;
  delete SizeXCmd;
  delete SizeYZCmd; 
  delete MagFieldCmd;
  delete TalNbCmd;
  delete TalDefCmd;
  delete TalPosiCmd;
  delete UpdateCmd;
  delete detDir;  
  delete testemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == MaterCmd )
   { Detector->SetMaterial(newValue);}
   
  if( command == SizeXCmd )
   { Detector->SetSizeX(SizeXCmd->GetNewDoubleValue(newValue));}
   
  if( command == SizeYZCmd )
   { Detector->SetSizeYZ(SizeYZCmd->GetNewDoubleValue(newValue));}
      
  if( command == MagFieldCmd )
   { Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
      
  if( command == TalNbCmd )
   { Detector->SetTallyNumber(TalNbCmd->GetNewIntValue(newValue));}
   
  if (command == TalDefCmd)
   {
     G4int num; G4double v1, v2, v3;
     G4String unt, mat;
     std::istringstream is(newValue);
     is >> num >> mat >> v1 >> v2 >> v3 >> unt;
     G4String material=mat;
     v1 *= G4UIcommand::ValueOf(unt);
     v2 *= G4UIcommand::ValueOf(unt);
     v3 *= G4UIcommand::ValueOf(unt);          
     Detector->SetTallyMaterial (num,material);
     Detector->SetTallySize(num,G4ThreeVector(v1,v2,v3));
   }
   
  if (command == TalPosiCmd)
   {
     G4int num; G4double v1, v2, v3;
     G4String unt;
     std::istringstream is(newValue);
     is >> num >> v1 >> v2 >> v3 >> unt;
     v1 *= G4UIcommand::ValueOf(unt);
     v2 *= G4UIcommand::ValueOf(unt);
     v3 *= G4UIcommand::ValueOf(unt);          
     Detector->SetTallyPosition(num,G4ThreeVector(v1,v2,v3));
   }      

  if( command == UpdateCmd )
   { Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
