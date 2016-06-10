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
// $Id: G4EnergyLossMessenger.cc 85424 2014-10-29 08:23:44Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4EnergyLossMessenger
//
// Author:        Michel Maire
//
// Creation date: 17-03-2011 (original version of 22-06-2000)
//
// Modifications:
// 10-01-06 SetStepLimits -> SetStepFunction (V.Ivanchenko)
// 10-01-06 PreciseRange -> CSDARange (V.Ivanchenko)
// 10-05-06 Add command MscStepLimit (V.Ivanchenko) 
// 10-10-06 Add DEDXBinning command (V.Ivanchenko)
// 07-02-07 Add MscLateralDisplacement command (V.Ivanchenko)
// 12-02-07 Add SetSkin, SetLinearLossLimit (V.Ivanchenko)
// 15-03-07 Send a message "/run/physicsModified" if reinitialisation
//          is needed after the command (V.Ivanchenko)
// 16-03-07 modify /process/eLoss/minsubsec command (V.Ivanchenko)
// 18-05-07 add /process/msc directory and commands (V.Ivanchenko)
// 11-03-08 add /process/em directory and commands (V.Ivanchenko)
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EnergyLossMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4EmProcessOptions.hh"
#include "G4UImanager.hh"
#include "G4MscStepLimitType.hh"
#include "G4EmProcessOptions.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossMessenger::G4EnergyLossMessenger()
{
  
  opt = 0;
  SubSecCmd = new G4UIcmdWithABool("/process/eLoss/subsec",this);
  SubSecCmd->SetGuidance("Switch true/false the subcutoff generation.");
  SubSecCmd->SetParameterName("choice",true);
  SubSecCmd->SetDefaultValue(true);
  SubSecCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  StepFuncCmd = new G4UIcommand("/process/eLoss/StepFunction",this);
  StepFuncCmd->SetGuidance("Set the energy loss step limitation parameters.");
  StepFuncCmd->SetGuidance("  dRoverR   : max Range variation per step");
  StepFuncCmd->SetGuidance("  finalRange: range for final step");

  G4UIparameter* dRoverRPrm = new G4UIparameter("dRoverR",'d',false);
  dRoverRPrm->SetGuidance("max Range variation per step (fractional number)");
  dRoverRPrm->SetParameterRange("dRoverR>0. && dRoverR<=1.");
  StepFuncCmd->SetParameter(dRoverRPrm);

  G4UIparameter* finalRangePrm = new G4UIparameter("finalRange",'d',false);
  finalRangePrm->SetGuidance("range for final step");
  finalRangePrm->SetParameterRange("finalRange>0.");
  StepFuncCmd->SetParameter(finalRangePrm);

  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',true);
  unitPrm->SetGuidance("unit of finalRange");
  unitPrm->SetDefaultValue("mm");
  G4String unitCandidates = 
    G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitCandidates);

  StepFuncCmd->SetParameter(unitPrm);
  StepFuncCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  IntegCmd = new G4UIcmdWithABool("/process/eLoss/integral",this);
  IntegCmd->SetGuidance("Switch true/false the integral option");
  IntegCmd->SetParameterName("integ",true);
  IntegCmd->SetDefaultValue(true);
  IntegCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  pixeXsCmd = new G4UIcmdWithAString("/process/em/pixeXSmodel",this);
  pixeXsCmd->SetGuidance("The name of PIXE cross section");
  pixeXsCmd->SetParameterName("pixeXS",true);
  pixeXsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  pixeeXsCmd = new G4UIcmdWithAString("/process/em/pixeElecXSmodel",this);
  pixeeXsCmd->SetGuidance("The name of PIXE cross section for electron");
  pixeeXsCmd->SetParameterName("pixeEXS",true);
  pixeeXsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  deexCmd = new G4UIcommand("/process/em/deexcitation",this);
  deexCmd->SetGuidance("Set deexcitation flags per G4Region.");
  deexCmd->SetGuidance("  regName   : G4Region name");
  deexCmd->SetGuidance("  flagFluo  : Fluorescence");
  deexCmd->SetGuidance("  flagAuger : Auger");
  deexCmd->SetGuidance("  flagPIXE  : PIXE");

  G4UIparameter* regName = new G4UIparameter("regName",'s',false);
  deexCmd->SetParameter(regName);

  G4UIparameter* flagFluo = new G4UIparameter("flagFluo",'s',false);
  deexCmd->SetParameter(flagFluo);

  G4UIparameter* flagAuger = new G4UIparameter("flagAuger",'s',false);
  deexCmd->SetParameter(flagAuger);

  G4UIparameter* flagPIXE = new G4UIparameter("flagPIXE",'s',false);
  deexCmd->SetParameter(flagPIXE);

  bfCmd = new G4UIcommand("/process/em/setBiasingFactor",this);
  bfCmd->SetGuidance("Set factor for the process cross section.");
  bfCmd->SetGuidance("  procName   : process name");
  bfCmd->SetGuidance("  procFact   : factor");
  bfCmd->SetGuidance("  flagFact   : flag to change weight");

  G4UIparameter* procName = new G4UIparameter("procName",'s',false);
  bfCmd->SetParameter(procName);

  G4UIparameter* procFact = new G4UIparameter("procFact",'d',false);
  bfCmd->SetParameter(procFact);

  G4UIparameter* flagFact = new G4UIparameter("flagFact",'s',false);
  bfCmd->SetParameter(flagFact);
  bfCmd->AvailableForStates(G4State_Idle);

  fiCmd = new G4UIcommand("/process/em/setForcedInteraction",this);
  fiCmd->SetGuidance("Set factor for the process cross section.");
  fiCmd->SetGuidance("  procNam    : process name");
  fiCmd->SetGuidance("  regNam     : region name");
  fiCmd->SetGuidance("  tlength    : fixed target length");
  fiCmd->SetGuidance("  tflag      : flag to change weight");

  G4UIparameter* procNam = new G4UIparameter("procNam",'s',false);
  fiCmd->SetParameter(procNam);

  G4UIparameter* regNam = new G4UIparameter("regNam",'s',false);
  fiCmd->SetParameter(regNam);

  G4UIparameter* tlength = new G4UIparameter("tlength",'d',false);
  fiCmd->SetParameter(tlength);

  G4UIparameter* unitT = new G4UIparameter("unitT",'s',true);
  fiCmd->SetParameter(unitT);
  unitT->SetGuidance("unit of tlength");

  G4UIparameter* flagT = new G4UIparameter("tflag",'s',true);
  fiCmd->SetParameter(flagT);
  fiCmd->AvailableForStates(G4State_Idle);

  brCmd = new G4UIcommand("/process/em/setSecBiasing",this);
  brCmd->SetGuidance("Set bremsstrahlung or delta-e- splitting/Russian roullette per region.");
  brCmd->SetGuidance("  bProcNam : process name");
  brCmd->SetGuidance("  bRegNam  : region name");
  brCmd->SetGuidance("  bFactor  : number of splitted gamma or probability of Russian roulette");
  brCmd->SetGuidance("  bEnergy  : max energy of a secondary for this biasing method");

  G4UIparameter* bProcNam = new G4UIparameter("bProcNam",'s',false);
  brCmd->SetParameter(bProcNam);

  G4UIparameter* bRegNam = new G4UIparameter("bRegNam",'s',false);
  brCmd->SetParameter(bRegNam);

  G4UIparameter* bFactor = new G4UIparameter("bFactor",'d',false);
  brCmd->SetParameter(bFactor);

  G4UIparameter* bEnergy = new G4UIparameter("bEnergy",'d',false);
  brCmd->SetParameter(bEnergy);

  G4UIparameter* bUnit = new G4UIparameter("bUnit",'s',true);
  brCmd->SetParameter(bUnit);
  brCmd->SetGuidance("unit of energy");

  brCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossMessenger::~G4EnergyLossMessenger()
{
  delete opt;
  delete SubSecCmd;
  delete StepFuncCmd;
  delete IntegCmd;
  delete pixeXsCmd;
  delete pixeeXsCmd;
  delete bfCmd;
  delete fiCmd;
  delete brCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(!opt) { opt = new G4EmProcessOptions(); }

  if(command == SubSecCmd) {
    opt->SetSubCutoff(SubSecCmd->GetNewBoolValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  } else if (command == StepFuncCmd) {
    G4double v1,v2;
    G4String unt;
    std::istringstream is(newValue);
    is >> v1 >> v2 >> unt;
    v2 *= G4UIcommand::ValueOf(unt);
    opt->SetStepFunction(v1,v2);
  } else if (command == deexCmd) {
    G4String s1 (""), s2(""), s3(""), s4("");
    G4bool b2(false), b3(false), b4(false);
    std::istringstream is(newValue);
    is >> s1 >> s2 >> s3 >> s4;
    if(s2 == "true") { b2 = true; }
    if(s3 == "true") { b3 = true; }
    if(s4 == "true") { b4 = true; }
    opt->SetDeexcitationActiveRegion(s1,b2,b3,b4);
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  } else if (command == pixeXsCmd) {
    G4String name;
    if (newValue == "ecpssr_analytical") 
      {name = "ECPSSR_Analytical";}
    else if (newValue == "ecpssr_interpolated") 
      {name = "ECPSSR_FormFactor";}
    else 
      {name = newValue;}
    opt->SetPIXECrossSectionModel(name);
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  } else if (command == pixeeXsCmd) {
    opt->SetPIXEElectronCrossSectionModel(newValue);
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  } else if (command == IntegCmd) {
    opt->SetIntegral(IntegCmd->GetNewBoolValue(newValue));
  } else if (command == bfCmd) {
    G4double v1(1.0);
    G4String s0(""),s1("");
    std::istringstream is(newValue);
    is >> s0 >> v1 >> s1;
    G4bool yes = false;
    if(s1 == "true") { yes = true; }
    opt->SetProcessBiasingFactor(s0,v1,yes);
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  } else if (command == fiCmd) {
    G4double v1(0.0);
    G4String s1(""),s2(""),s3(""),unt("mm");
    std::istringstream is(newValue);
    is >> s1 >> s2 >> v1 >> unt >> s3;
    G4bool yes = false;
    if(s3 == "true") { yes = true; }
    v1 *= G4UIcommand::ValueOf(unt);
    opt->ActivateForcedInteraction(s1,v1,s2,yes);
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  } else if (command == brCmd) {
    G4double fb(1.0),en(1.e+30);
    G4String s1(""),s2(""),unt("MeV");
    std::istringstream is(newValue);
    is >> s1 >> s2 >> fb >> en >> unt;
    en *= G4UIcommand::ValueOf(unt);    
    if (s1=="phot"||s1=="compt"||s1=="conv") 
                opt->ActivateSecondaryBiasingForGamma(s1,s2,fb,en);
    else opt->ActivateSecondaryBiasing(s1,s2,fb,en);
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
