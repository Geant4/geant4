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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4EmLowEParametersMessenger
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 07-05-2019 
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmLowEParametersMessenger.hh"
#include "G4EmFluoDirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImanager.hh"
#include "G4EmLowEParameters.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmLowEParametersMessenger::G4EmLowEParametersMessenger(G4EmLowEParameters* ptr) 
  : theParameters(ptr)
{
  deCmd = new G4UIcmdWithABool("/process/em/fluo",this);
  deCmd->SetGuidance("Enable/disable atomic deexcitation");
  deCmd->SetParameterName("fluoFlag",true);
  deCmd->SetDefaultValue(false);
  deCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);
  deCmd->SetToBeBroadcasted(false);

  dirFluoCmd = new G4UIcmdWithABool("/process/em/fluoBearden",this);
  dirFluoCmd->SetGuidance("Enable/disable usage of Bearden fluorescence files");
  dirFluoCmd->SetParameterName("fluoBeardenFlag",true);
  dirFluoCmd->SetDefaultValue(false);
  dirFluoCmd->AvailableForStates(G4State_PreInit,G4State_Init);
  dirFluoCmd->SetToBeBroadcasted(false);

  dirFluoCmd1 = new G4UIcmdWithABool("/process/em/fluoANSTO",this);
  dirFluoCmd1->SetGuidance("Enable/disable usage of ANSTO fluorescence files");
  dirFluoCmd1->SetParameterName("fluoANSTOFlag",true);
  dirFluoCmd1->SetDefaultValue(false);
  dirFluoCmd1->AvailableForStates(G4State_PreInit,G4State_Init);
  dirFluoCmd1->SetToBeBroadcasted(false);

  auCmd = new G4UIcmdWithABool("/process/em/auger",this);
  auCmd->SetGuidance("Enable/disable Auger electrons production");
  auCmd->SetParameterName("augerFlag",true);
  auCmd->SetDefaultValue(false);
  auCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);
  auCmd->SetToBeBroadcasted(false);

  auCascadeCmd = new G4UIcmdWithABool("/process/em/augerCascade",this);
  auCascadeCmd->SetGuidance("Enable/disable simulation of cascade of Auger electrons");
  auCascadeCmd->SetParameterName("augerCascadeFlag",true);
  auCascadeCmd->SetDefaultValue(false);
  auCascadeCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);
  auCascadeCmd->SetToBeBroadcasted(false);

  pixeCmd = new G4UIcmdWithABool("/process/em/pixe",this);
  pixeCmd->SetGuidance("Enable/disable PIXE simulation");
  pixeCmd->SetParameterName("pixeFlag",true);
  pixeCmd->SetDefaultValue(false);
  pixeCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);
  pixeCmd->SetToBeBroadcasted(false);

  dcutCmd = new G4UIcmdWithABool("/process/em/deexcitationIgnoreCut",this);
  dcutCmd->SetGuidance("Enable/Disable usage of cuts in de-excitation module");
  dcutCmd->SetParameterName("deexcut",true);
  dcutCmd->SetDefaultValue(false);
  dcutCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);
  dcutCmd->SetToBeBroadcasted(false);

  dnafCmd = new G4UIcmdWithABool("/process/dna/UseDNAFast",this);
  dnafCmd->SetGuidance("Enable usage of fast sampling for DNA models");
  dnafCmd->SetParameterName("dnaf",true);
  dnafCmd->SetDefaultValue(false);
  dnafCmd->AvailableForStates(G4State_PreInit);
  dnafCmd->SetToBeBroadcasted(false);

  dnasCmd = new G4UIcmdWithABool("/process/dna/UseDNAStationary",this);
  dnasCmd->SetGuidance("Enable usage of Stationary option for DNA models");
  dnasCmd->SetParameterName("dnas",true);
  dnasCmd->SetDefaultValue(false);
  dnasCmd->AvailableForStates(G4State_PreInit);
  dnasCmd->SetToBeBroadcasted(false);

  dnamscCmd = new G4UIcmdWithABool("/process/dna/UseDNAElectronMsc",this);
  dnamscCmd->SetGuidance("Enable usage of e- msc for DNA");
  dnamscCmd->SetParameterName("dnamsc",true);
  dnamscCmd->SetDefaultValue(false);
  dnamscCmd->AvailableForStates(G4State_PreInit);
  dnamscCmd->SetToBeBroadcasted(false);

  direFluoCmd = new G4UIcmdWithAString("/process/em/fluoDirectory",this);
  direFluoCmd->SetGuidance("The name of PIXE cross section");
  direFluoCmd->SetParameterName("fluoDirectory",true);
  direFluoCmd->SetCandidates("Default Bearden ANSTO XDB_EADL");
  direFluoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  direFluoCmd->SetToBeBroadcasted(false);

  pixeXsCmd = new G4UIcmdWithAString("/process/em/pixeXSmodel",this);
  pixeXsCmd->SetGuidance("The name of PIXE cross section");
  pixeXsCmd->SetParameterName("pixeXS",true);
  pixeXsCmd->SetCandidates("ECPSSR_Analytical Empirical ECPSSR_FormFactor ECPSSR_ANSTO");
  pixeXsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  pixeXsCmd->SetToBeBroadcasted(false);

  pixeeXsCmd = new G4UIcmdWithAString("/process/em/pixeElecXSmodel",this);
  pixeeXsCmd->SetGuidance("The name of PIXE cross section for electron");
  pixeeXsCmd->SetParameterName("pixeEXS",true);
  pixeeXsCmd->SetCandidates("ECPSSR_Analytical Empirical Livermore Penelope");
  pixeeXsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  pixeeXsCmd->SetToBeBroadcasted(false);

  livCmd = new G4UIcmdWithAString("/process/em/LivermoreData",this);
  livCmd->SetGuidance("The name of Livermore data directory");
  livCmd->SetParameterName("livDir",true);
  livCmd->SetCandidates("livermore epics_2017");
  livCmd->AvailableForStates(G4State_PreInit);
  livCmd->SetToBeBroadcasted(false);

  dnaSolCmd = new G4UIcmdWithAString("/process/dna/e-SolvationSubType",this);
  dnaSolCmd->SetGuidance("The name of e- solvation DNA model");
  dnaSolCmd->SetParameterName("dnaSol",true);
  dnaSolCmd->SetCandidates("Ritchie1994 Terrisol1990 Meesungnoen2002 Kreipl2009 Meesungnoen2002_amorphous");
  dnaSolCmd->AvailableForStates(G4State_PreInit);
  dnaSolCmd->SetToBeBroadcasted(false);

  meCmd = new G4UIcmdWithAString("/process/em/AddMicroElecRegion",this);
  meCmd->SetGuidance("Activate MicroElec model in the G4Region");
  meCmd->SetParameterName("MicroElec",true);
  meCmd->AvailableForStates(G4State_PreInit);
  meCmd->SetToBeBroadcasted(false);

  dnaCmd = new G4UIcommand("/process/em/AddDNARegion",this);
  dnaCmd->SetGuidance("Activate DNA in a G4Region.");
  dnaCmd->SetGuidance("  regName   : G4Region name");
  dnaCmd->SetGuidance("  dnaType   : DNA_opt0, DNA_Opt2, DNA_Opt4, DNA_Opt4a, DNA_Opt6, DNA_Opt6a, DNA_Opt7");
  dnaCmd->AvailableForStates(G4State_PreInit);
  dnaCmd->SetToBeBroadcasted(false);

  auto regName = new G4UIparameter("regName",'s',false);
  dnaCmd->SetParameter(regName);

  auto type = new G4UIparameter("dnaType",'s',false);
  dnaCmd->SetParameter(type);
  type->SetParameterCandidates("DNA_Opt0 DNA_Opt2 DNA_Opt4 DNA_Opt4a DNA_Opt6 DNA_Opt6a DNA_Opt7");

  deexCmd = new G4UIcommand("/process/em/deexcitation",this);
  deexCmd->SetGuidance("Set deexcitation flags per G4Region.");
  deexCmd->SetGuidance("  regName   : G4Region name");
  deexCmd->SetGuidance("  flagFluo  : Fluorescence");
  deexCmd->SetGuidance("  flagAuger : Auger");
  deexCmd->SetGuidance("  flagPIXE  : PIXE");
  deexCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);
  deexCmd->SetToBeBroadcasted(false);

  auto regNameD = new G4UIparameter("regName",'s',false);
  deexCmd->SetParameter(regNameD);

  auto flagFluo = new G4UIparameter("flagFluo",'s',false);
  deexCmd->SetParameter(flagFluo);

  auto flagAuger = new G4UIparameter("flagAuger",'s',false);
  deexCmd->SetParameter(flagAuger);

  auto flagPIXE = new G4UIparameter("flagPIXE",'s',false);
  deexCmd->SetParameter(flagPIXE);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmLowEParametersMessenger::~G4EmLowEParametersMessenger()
{
  delete deCmd;
  delete dirFluoCmd;
  delete dirFluoCmd1;
  delete auCmd;
  delete auCascadeCmd;
  delete pixeCmd;
  delete dcutCmd;
  delete dnafCmd;
  delete dnasCmd;
  delete dnamscCmd;
  delete pixeXsCmd;
  delete pixeeXsCmd;
  delete livCmd;
  delete dnaSolCmd;
  delete direFluoCmd;
  delete meCmd;
  delete dnaCmd;
  delete deexCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmLowEParametersMessenger::SetNewValue(G4UIcommand* command, 
                                              G4String newValue)
{
  G4bool physicsModified = false;
  if (command == deCmd) {
    theParameters->SetFluo(deCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == dirFluoCmd) {
    theParameters->SetBeardenFluoDir(dirFluoCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == dirFluoCmd1) {
    theParameters->SetANSTOFluoDir(dirFluoCmd1->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == auCmd) {
    theParameters->SetAuger(auCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == auCascadeCmd) {
    theParameters->SetAuger(auCascadeCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == pixeCmd) {
    theParameters->SetPixe(pixeCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == dcutCmd) {
    theParameters->SetDeexcitationIgnoreCut(dcutCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == dnafCmd) {
    theParameters->SetDNAFast(dnafCmd->GetNewBoolValue(newValue));
  } else if (command == dnasCmd) {
    theParameters->SetDNAStationary(dnasCmd->GetNewBoolValue(newValue));
  } else if (command == dnamscCmd) {
    theParameters->SetDNAElectronMsc(dnamscCmd->GetNewBoolValue(newValue));
  } else if (command == dnaSolCmd) {
    G4DNAModelSubType ttt = fDNAUnknownModel;
    if(newValue == "Ritchie1994") { 
      ttt = fRitchie1994eSolvation; 
    } else if(newValue == "Terrisol1990") { 
      ttt = fTerrisol1990eSolvation; 
    } else if (newValue == "Meesungnoen2002") { 
      ttt = fMeesungnoen2002eSolvation;
    } else if (newValue == "Meesungnoen2002_amorphous") {
      ttt = fMeesungnoensolid2002eSolvation;
    } else if (newValue == "Kreipl2009") {
      ttt = fKreipl2009eSolvation;
    }
    theParameters->SetDNAeSolvationSubType(ttt);
  } else if (command == direFluoCmd) {
    G4EmFluoDirectory ttt = fluoDefault;
    if(newValue == "Bearden") { ttt = fluoBearden; }
    else if(newValue == "ANSTO") { ttt = fluoANSTO; }
    else if(newValue == "XDB_EADL") { ttt = fluoXDB_EADL; }
    theParameters->SetFluoDirectory(ttt);
  } else if (command == pixeXsCmd) {
    theParameters->SetPIXECrossSectionModel(newValue);
    physicsModified = true;
  } else if (command == pixeeXsCmd) {
    theParameters->SetPIXEElectronCrossSectionModel(newValue);
    physicsModified = true;
  } else if (command == livCmd) {
    theParameters->SetLivermoreDataDir(newValue);
  } else if (command == meCmd) {
    theParameters->AddMicroElec(newValue);
  } else if (command == dnaCmd) {
    G4String s1(""),s2("");
    std::istringstream is(newValue);
    is >> s1 >> s2;
    theParameters->AddDNA(s1, s2);
  } else if (command == deexCmd) {
    G4String s1 (""), s2(""), s3(""), s4("");
    G4bool b2(false), b3(false), b4(false);
    std::istringstream is(newValue);
    is >> s1 >> s2 >> s3 >> s4;
    if(s2 == "true") { b2 = true; }
    if(s3 == "true") { b3 = true; }
    if(s4 == "true") { b4 = true; }
    theParameters->SetDeexActiveRegion(s1,b2,b3,b4);
    physicsModified = true;
  }
  
  if(physicsModified) {
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
