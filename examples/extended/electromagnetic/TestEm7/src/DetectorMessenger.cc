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
/// \file electromagnetic/TestEm7/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
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
:G4UImessenger(),fDetector(Det),
 fTestemDir(nullptr),
 fDetDir(nullptr),    
 fMaterCmd(nullptr),
 fWMaterCmd(nullptr),
 fSizeXCmd(nullptr),
 fSizeYZCmd(nullptr),    
 fMagFieldCmd(nullptr),
 fTalNbCmd(nullptr),    
 fTalDefCmd(nullptr),
 fTalPosiCmd(nullptr)
{ 
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance(" detector control.");
  
  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction commands");
      
  fMaterCmd = new G4UIcmdWithAString("/testem/det/setMat",this);
  fMaterCmd->SetGuidance("Select material of the box.");
  fMaterCmd->SetParameterName("choice",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fWMaterCmd = new G4UIcmdWithAString("/testem/det/setWorldMat",this);
  fWMaterCmd->SetGuidance("Select material of the world.");
  fWMaterCmd->SetParameterName("choice",false);
  fWMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSizeXCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeX",this);
  fSizeXCmd->SetGuidance("Set sizeX of the absorber");
  fSizeXCmd->SetParameterName("SizeX",false);
  fSizeXCmd->SetRange("SizeX>0.");
  fSizeXCmd->SetUnitCategory("Length");
  fSizeXCmd->AvailableForStates(G4State_PreInit);
  
  fSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeYZ",this);
  fSizeYZCmd->SetGuidance("Set sizeYZ of the absorber");
  fSizeYZCmd->SetParameterName("SizeYZ",false);
  fSizeYZCmd->SetRange("SizeYZ>0.");
  fSizeYZCmd->SetUnitCategory("Length");
  fSizeYZCmd->AvailableForStates(G4State_PreInit);
        
  fMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);  
  fMagFieldCmd->SetGuidance("Define magnetic field.");
  fMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMagFieldCmd->SetParameterName("Bz",false);
  fMagFieldCmd->SetUnitCategory("Magnetic flux density");
  fMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fTalNbCmd = new G4UIcmdWithAnInteger("/testem/det/tallyNumber",this);
  fTalNbCmd->SetGuidance("Set number of fTallies.");
  fTalNbCmd->SetParameterName("tallyNb",false);
  fTalNbCmd->SetRange("tallyNb>=0");
  fTalNbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fTalDefCmd = new G4UIcommand("/testem/det/tallyDefinition",this);
  fTalDefCmd->SetGuidance("Set tally nb, box dimensions.");
  fTalDefCmd->SetGuidance("  tally number : from 0 to tallyNumber");
  fTalDefCmd->SetGuidance("  material name");
  fTalDefCmd->SetGuidance("  dimensions (3-vector with unit)");
  //
  G4UIparameter* fTalNbPrm = new G4UIparameter("tallyNb",'i',false);
  fTalNbPrm->SetGuidance("tally number : from 0 to tallyNumber");
  fTalNbPrm->SetParameterRange("tallyNb>=0");
  fTalDefCmd->SetParameter(fTalNbPrm);
  //
  G4UIparameter* SizeXPrm = new G4UIparameter("sizeX",'d',false);
  SizeXPrm->SetGuidance("sizeX");
  SizeXPrm->SetParameterRange("sizeX>0.");
  fTalDefCmd->SetParameter(SizeXPrm);
  //    
  G4UIparameter* SizeYPrm = new G4UIparameter("sizeY",'d',false);
  SizeYPrm->SetGuidance("sizeY");
  SizeYPrm->SetParameterRange("sizeY>0.");
  fTalDefCmd->SetParameter(SizeYPrm);
  //    
  G4UIparameter* SizeZPrm = new G4UIparameter("sizeZ",'d',false);
  SizeZPrm->SetGuidance("sizeZ");
  SizeZPrm->SetParameterRange("sizeZ>0.");
  fTalDefCmd->SetParameter(SizeZPrm);    
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of dimensions");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitList);
  fTalDefCmd->SetParameter(unitPrm);
  //
  fTalDefCmd->AvailableForStates(G4State_PreInit);

  fTalPosiCmd = new G4UIcommand("/testem/det/tallyPosition",this);
  fTalPosiCmd->SetGuidance("Set tally nb, position");
  fTalPosiCmd->SetGuidance("  tally number : from 0 to tallyNumber");
  fTalPosiCmd->SetGuidance("  position (3-vector with unit)");
  //
  G4UIparameter* fTalNumPrm = new G4UIparameter("tallyNum",'i',false);
  fTalNumPrm->SetGuidance("tally number : from 0 to tallyNumber");
  fTalNumPrm->SetParameterRange("tallyNum>=0");
  fTalPosiCmd->SetParameter(fTalNumPrm);
  //    
  G4UIparameter* PosiXPrm = new G4UIparameter("posiX",'d',false);
  PosiXPrm->SetGuidance("position X");
  fTalPosiCmd->SetParameter(PosiXPrm);
  //    
  G4UIparameter* PosiYPrm = new G4UIparameter("posiY",'d',false);
  PosiYPrm->SetGuidance("position Y");
  fTalPosiCmd->SetParameter(PosiYPrm);
  //
  G4UIparameter* PosiZPrm = new G4UIparameter("posiZ",'d',false);
  PosiZPrm->SetGuidance("position Z");
  fTalPosiCmd->SetParameter(PosiZPrm);      
  //
  G4UIparameter* unitPr = new G4UIparameter("unit",'s',false);
  unitPr->SetGuidance("unit of position");
  unitPr->SetParameterCandidates(unitList);
  fTalPosiCmd->SetParameter(unitPr);
  //
  fTalPosiCmd->AvailableForStates(G4State_PreInit);
        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterCmd;
  delete fWMaterCmd;
  delete fSizeXCmd;
  delete fSizeYZCmd; 
  delete fMagFieldCmd;
  delete fTalNbCmd;
  delete fTalDefCmd;
  delete fTalPosiCmd;
  delete fDetDir;  
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fMaterCmd )
   { fDetector->SetMaterial(newValue);}
   
  if( command == fWMaterCmd )
   { fDetector->SetWorldMaterial(newValue);}
   
  if( command == fSizeXCmd )
   { fDetector->SetSizeX(fSizeXCmd->GetNewDoubleValue(newValue));}
   
  if( command == fSizeYZCmd )
   { fDetector->SetSizeYZ(fSizeYZCmd->GetNewDoubleValue(newValue));}
      
  if( command == fMagFieldCmd )
   { fDetector->SetMagField(fMagFieldCmd->GetNewDoubleValue(newValue));}
      
  if( command == fTalNbCmd )
   { fDetector->SetTallyNumber(fTalNbCmd->GetNewIntValue(newValue));}
   
  if (command == fTalDefCmd)
   {
     G4int num; G4double v1, v2, v3;
     G4String unt;
     std::istringstream is(newValue);
     is >> num >> v1 >> v2 >> v3 >> unt;
     G4ThreeVector vec(v1,v2,v3);
     vec *= G4UIcommand::ValueOf(unt);
     fDetector->SetTallySize(num,vec);
   }
   
  if (command == fTalPosiCmd)
   {
     G4int num; G4double v1, v2, v3;
     G4String unt;
     std::istringstream is(newValue);
     is >> num >> v1 >> v2 >> v3 >> unt;
     G4ThreeVector vec(v1,v2,v3);
     vec *= G4UIcommand::ValueOf(unt);
     fDetector->SetTallyPosition(num,vec);
   }      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
