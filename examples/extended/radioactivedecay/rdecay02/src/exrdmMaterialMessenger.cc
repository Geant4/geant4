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
/// \file radioactivedecay/rdecay02/src/exrdmMaterialMessenger.cc
/// \brief Implementation of the exrdmMaterialMessenger class
//
////////////////////////////////////////////////////////////////////////////////
//
#include "exrdmMaterialMessenger.hh"
#include "exrdmMaterial.hh"

#include "G4SystemOfUnits.hh"

#include <sstream>

////////////////////////////////////////////////////////////////////////////////
//
exrdmMaterialMessenger::exrdmMaterialMessenger (exrdmMaterial * exrdmMat)
  :fMaterialsManager(exrdmMat)
{ 
  fMaterialDir = new G4UIdirectory("/geometry/material/");
  fMaterialDir->SetGuidance(" Controls for defining geometry materials" );

  fAddCmd = new G4UIcommand("/geometry/material/add",this);
  fAddCmd->SetGuidance(
    "  add a mateial by name, composition formula and density");
  fAddCmd->SetGuidance("  name: e.g. water ");
  fAddCmd->SetGuidance("  formula (e.g. H2-O for water");
  fAddCmd->SetGuidance("  density (in units of g/cm3) : den>0.");
  G4UIparameter* MatName = new G4UIparameter("material",'s',false);
  MatName->SetGuidance("material name");
  fAddCmd->SetParameter(MatName);
  //
  G4UIparameter* MatForm = new G4UIparameter("formula",'s',false);
  MatForm->SetGuidance("material formula");
  fAddCmd->SetParameter(MatForm);
  //    
  G4UIparameter* DenPrm = new G4UIparameter("density",'d',false);
  DenPrm->SetGuidance("density of the material");
  DenPrm->SetParameterRange("density >0.");
  fAddCmd->SetParameter(DenPrm);
  fAddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* StatePrm = new G4UIparameter("state",'s',true);
  StatePrm->SetGuidance("state of the material (optional): gas | solid");
  fAddCmd->SetParameter(StatePrm);
  fAddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* TempPrm = new G4UIparameter("temp",'d',true);
  TempPrm->SetGuidance("temperature of the material in Kelvin (optional)");
  fAddCmd->SetParameter(TempPrm);
  fAddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* PresPrm = new G4UIparameter("pres",'d',true);
  PresPrm->SetGuidance("pressure of the gas material in Pascal (optional)");
  fAddCmd->SetParameter(PresPrm);
  fAddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //
  fDeleteIntCmd = new G4UIcmdWithAnInteger("/geometry/material/delete",this);
  fDeleteIntCmd->SetGuidance("Delete material by its index");
  fDeleteIntCmd->SetParameterName("matIdx",false);
  fDeleteIntCmd->SetRange("matIdx>=0 && matIdx<100");
  fDeleteIntCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fDeleteNameCmd = new G4UIcmdWithAString("/geometry/material/deleteName",this);
  fDeleteNameCmd->SetGuidance("Delete material by its name.");
  fDeleteNameCmd->SetParameterName("DeleteName",false);
  fDeleteNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  fListCmd = new G4UIcmdWithoutParameter("/geometry/material/list",this);
  fListCmd->SetGuidance("List the materials defined");
  fListCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}
////////////////////////////////////////////////////////////////////////////////
//
exrdmMaterialMessenger::~exrdmMaterialMessenger ()
{
  delete fMaterialDir;
  delete fAddCmd;
  delete fDeleteIntCmd;
  delete fDeleteNameCmd;
  delete fListCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
void exrdmMaterialMessenger::SetNewValue (G4UIcommand* command,
                                          G4String newValue)
{    
  if (command == fDeleteIntCmd) {
    fMaterialsManager->DeleteMaterial(fDeleteIntCmd->GetNewIntValue(newValue));

  } else if (command == fDeleteNameCmd) {
    fMaterialsManager->DeleteMaterial(newValue);

  } else if (command == fListCmd) {
    fMaterialsManager->ListMaterial();

  } else if (command == fAddCmd) {
    G4double den, tem, pres ;
    G4String state;
    char mat[80], form[80], stat[10];
    stat[0] = ' ';
    tem = pres = -1.;
    const char* t = newValue;
    std::istringstream is(t);
    is >>mat >>form >>den >>stat >> tem >> pres ;
    G4String material=mat;
    G4String formula=form;
    if (pres == -1.) { 
      state = "";
    } else {
      state = stat;
    }
    //    G4cout<< "stat = " <<state<< "tem = " << tem<< " pre = " << pres << G4endl;
    //     tick *= G4UIcommand::ValueOf(unt);
    fMaterialsManager->AddMaterial(material,formula,den*g/cm3,state,tem,pres);
  }
}
////////////////////////////////////////////////////////////////////////////////
