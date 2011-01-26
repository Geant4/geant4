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
//
// $Id: G4EnergyLossMessenger.cc,v 1.41 2011-01-03 19:34:03 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4EnergyLossMessenger
//
// Author:        Michel Maire
//
// Creation date: 22-06-2000
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

#include "G4VEnergyLoss.hh"

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
  eLossDirectory = new G4UIdirectory("/process/eLoss/");
  eLossDirectory->SetGuidance("Commands for EM processes.");
  mscDirectory = new G4UIdirectory("/process/msc/");
  mscDirectory->SetGuidance("Commands for EM scattering processes.");
  emDirectory = new G4UIdirectory("/process/em/");
  emDirectory->SetGuidance("General commands for EM processes.");

  RndmStepCmd = new G4UIcmdWithABool("/process/eLoss/rndmStep",this);
  RndmStepCmd->SetGuidance("Randomize the proposed step by eLoss.");
  RndmStepCmd->SetParameterName("choice",true);
  RndmStepCmd->SetDefaultValue(false);
  RndmStepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  EnlossFlucCmd = new G4UIcmdWithABool("/process/eLoss/fluct",this);
  EnlossFlucCmd->SetGuidance("Switch true/false the energy loss fluctuations.");
  EnlossFlucCmd->SetParameterName("choice",true);
  EnlossFlucCmd->SetDefaultValue(true);
  EnlossFlucCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SubSecCmd = new G4UIcmdWithABool("/process/eLoss/subsec",this);
  SubSecCmd->SetGuidance("Switch true/false the subcutoff generation.");
  SubSecCmd->SetParameterName("choice",true);
  SubSecCmd->SetDefaultValue(true);
  SubSecCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MinSubSecCmd = new G4UIcmdWithADouble("/process/eLoss/minsubsec",this);
  MinSubSecCmd->SetGuidance("Set the ratio subcut/cut ");
  MinSubSecCmd->SetParameterName("rcmin",true);
  MinSubSecCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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

  MinEnCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/minKinEnergy",this);
  MinEnCmd->SetGuidance("Set the min kinetic energy");
  MinEnCmd->SetParameterName("emin",true);
  MinEnCmd->SetUnitCategory("Energy");
  MinEnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MaxEnCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/maxKinEnergy",this);
  MaxEnCmd->SetGuidance("Set the max kinetic energy");
  MaxEnCmd->SetParameterName("emax",true);
  MaxEnCmd->SetUnitCategory("Energy");
  MaxEnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  IntegCmd = new G4UIcmdWithABool("/process/eLoss/integral",this);
  IntegCmd->SetGuidance("Switch true/false the integral option");
  IntegCmd->SetParameterName("integ",true);
  IntegCmd->SetDefaultValue(true);
  IntegCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  rangeCmd = new G4UIcmdWithABool("/process/eLoss/CSDARange",this);
  rangeCmd->SetGuidance("Switch true/false the CSDA range calculation");
  rangeCmd->SetParameterName("range",true);
  rangeCmd->SetDefaultValue(true);
  rangeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  lpmCmd = new G4UIcmdWithABool("/process/eLoss/LPM",this);
  lpmCmd->SetGuidance("The flag of the LPM effect calculation");
  lpmCmd->SetParameterName("lpm",true);
  lpmCmd->SetDefaultValue(true);
  lpmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  splCmd = new G4UIcmdWithABool("/process/em/spline",this);
  splCmd->SetGuidance("The flag of usage spline for Physics Vectors");
  splCmd->SetParameterName("spl",true);
  splCmd->SetDefaultValue(false);
  splCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  aplCmd = new G4UIcmdWithABool("/process/em/applyCuts",this);
  aplCmd->SetGuidance("The flag to Apply Cuts for gamma processes");
  aplCmd->SetParameterName("apl",true);
  aplCmd->SetDefaultValue(false);
  aplCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  deCmd = new G4UIcmdWithABool("/process/em/fluo",this);
  deCmd->SetGuidance("The flag to enable/disable deexcitation");
  deCmd->SetParameterName("fluoFlag",true);
  deCmd->SetDefaultValue(false);
  deCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  auCmd = new G4UIcmdWithABool("/process/em/auger",this);
  auCmd->SetGuidance("The flag to enable/disable Auger electrons");
  auCmd->SetParameterName("augerFlag",true);
  auCmd->SetDefaultValue(false);
  auCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  pixeCmd = new G4UIcmdWithABool("/process/em/pixe",this);
  pixeCmd->SetGuidance("The flag to enable/disable PIXE");
  pixeCmd->SetParameterName("pixeFlag",true);
  pixeCmd->SetDefaultValue(false);
  pixeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  pixeXsCmd = new G4UIcmdWithAString("/process/em/pixeXSmodel",this);
  pixeXsCmd->SetGuidance("The name of PIXE cross section");
  pixeXsCmd->SetParameterName("pixeXS",true);
  pixeXsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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

  dedxCmd = new G4UIcmdWithAnInteger("/process/eLoss/binsDEDX",this);
  dedxCmd->SetGuidance("Set number of bins for DEDX tables");
  dedxCmd->SetParameterName("binsDEDX",true);
  dedxCmd->SetDefaultValue(77);
  dedxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  lamCmd = new G4UIcmdWithAnInteger("/process/eLoss/binsLambda",this);
  lamCmd->SetGuidance("Set number of bins for Lambda tables");
  lamCmd->SetParameterName("binsL",true);
  lamCmd->SetDefaultValue(77);
  lamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  verCmd = new G4UIcmdWithAnInteger("/process/eLoss/verbose",this);
  verCmd->SetGuidance("Set verbose level for EM physics");
  verCmd->SetParameterName("verb",true);
  verCmd->SetDefaultValue(1);
  verCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ver1Cmd = new G4UIcmdWithAnInteger("/process/em/verbose",this);
  ver1Cmd->SetGuidance("Set verbose level for EM physics");
  ver1Cmd->SetParameterName("verb1",true);
  ver1Cmd->SetDefaultValue(1);
  ver1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  lllCmd = new G4UIcmdWithADouble("/process/eLoss/linLossLimit",this);
  lllCmd->SetGuidance("Set linearLossLimit parameter");
  lllCmd->SetParameterName("linlim",true);
  lllCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  labCmd = new G4UIcmdWithADouble("/process/eLoss/LambdaFactor",this);
  labCmd->SetGuidance("Set lambdaFactor parameter for integral option");
  labCmd->SetParameterName("Fl",true);
  labCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  mscCmd = new G4UIcmdWithAString("/process/msc/StepLimit",this);
  mscCmd->SetGuidance("Set msc step limitation type");
  mscCmd->SetParameterName("StepLim",true);
  mscCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  latCmd = new G4UIcmdWithABool("/process/msc/LateralDisplacement",this);
  latCmd->SetGuidance("Set flag of sampling of lateral displacement");
  latCmd->SetParameterName("lat",true);
  latCmd->SetDefaultValue(true);
  latCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  frCmd = new G4UIcmdWithADouble("/process/msc/RangeFactor",this);
  frCmd->SetGuidance("Set RangeFactor parameter for msc processes");
  frCmd->SetParameterName("Fr",true);
  frCmd->SetRange("Fr>0");
  frCmd->SetDefaultValue(0.04);
  frCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fgCmd = new G4UIcmdWithADouble("/process/msc/GeomFactor",this);
  fgCmd->SetGuidance("Set GeomFactor parameter for msc processes");
  fgCmd->SetParameterName("Fg",true);
  fgCmd->SetRange("Fg>0");
  fgCmd->SetDefaultValue(3.5);
  fgCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  mscfCmd = new G4UIcmdWithADouble("/process/msc/FactorForAngleLimit",this);
  mscfCmd->SetGuidance("Set factor for computation of a limit for -t (invariant trasfer)");
  mscfCmd->SetParameterName("Fact",true);
  mscfCmd->SetRange("Fact>0");
  mscfCmd->SetDefaultValue(1.);
  mscfCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  skinCmd = new G4UIcmdWithADouble("/process/msc/Skin",this);
  skinCmd->SetGuidance("Set skin parameter for msc processes");
  skinCmd->SetParameterName("skin",true);
  skinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  angCmd = new G4UIcmdWithADoubleAndUnit("/process/msc/ThetaLimit",this);
  angCmd->SetGuidance("Set the limit on the polar angle for msc and single scattering");
  angCmd->SetParameterName("theta",true);
  angCmd->SetUnitCategory("Angle");
  angCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossMessenger::~G4EnergyLossMessenger()
{
  delete opt;
  delete RndmStepCmd;
  delete EnlossFlucCmd;
  delete SubSecCmd;
  delete MinSubSecCmd;
  delete StepFuncCmd;
  delete deexCmd;
  delete eLossDirectory;
  delete mscDirectory;
  delete emDirectory;
  delete MinEnCmd;
  delete MaxEnCmd;
  delete IntegCmd;
  delete rangeCmd;
  delete lpmCmd;
  delete splCmd;
  delete aplCmd;
  delete latCmd;
  delete verCmd;
  delete ver1Cmd;
  delete mscCmd;
  delete dedxCmd;
  delete deCmd;
  delete auCmd;
  delete pixeCmd;
  delete pixeXsCmd;
  delete frCmd;
  delete fgCmd;
  delete lllCmd;
  delete lamCmd;
  delete labCmd;
  delete skinCmd;
  delete angCmd;
  delete mscfCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(!opt) opt = new G4EmProcessOptions();

  if (command == RndmStepCmd) {
    G4VEnergyLoss::SetRndmStep(RndmStepCmd->GetNewBoolValue(newValue));
    opt->SetRandomStep(RndmStepCmd->GetNewBoolValue(newValue));
    return;
  }
  if (command == EnlossFlucCmd) {
    G4VEnergyLoss::SetEnlossFluc(EnlossFlucCmd->GetNewBoolValue(newValue));
    opt->SetLossFluctuations(EnlossFlucCmd->GetNewBoolValue(newValue));
    return;
  }
  if (command == SubSecCmd) {
    G4VEnergyLoss::SetSubSec(SubSecCmd->GetNewBoolValue(newValue));
    opt->SetSubCutoff(SubSecCmd->GetNewBoolValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == MinSubSecCmd) {
    opt->SetMinSubRange(MinSubSecCmd->GetNewDoubleValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == StepFuncCmd) {
    G4double v1,v2;
    G4String unt;
    std::istringstream is(newValue);
    is >> v1 >> v2 >> unt;
    v2 *= G4UIcommand::ValueOf(unt);
    G4VEnergyLoss::SetStepFunction(v1,v2);
    opt->SetStepFunction(v1,v2);
    return;
  }
  if (command == deexCmd) {
    G4String s1 (""), s2(""), s3(""), s4("");
    G4bool b2(false), b3(false), b4(false);
    std::istringstream is(newValue);
    is >> s1 >> s2 >> s3 >> s4;
    if(s2 == "true") { b2 = true; }
    if(s3 == "true") { b3 = true; }
    if(s4 == "true") { b4 = true; }
    opt->SetDeexcitationActiveRegion(s1,b2,b3,b4);
    return;
  }
  if (command == deCmd) {
    opt->SetDeexcitationActive(deCmd->GetNewBoolValue(newValue));
    return;
  }
  if (command == auCmd) {
    opt->SetAugerActive(auCmd->GetNewBoolValue(newValue));
    return;
  }
  if (command == pixeCmd) {
    opt->SetPIXEActive(pixeCmd->GetNewBoolValue(newValue));
    return;
  }
  if (command == pixeXsCmd) {
    opt->SetPIXECrossSectionModel(newValue);
    return;
  }
  if (command == mscCmd) {
    if(newValue == "Minimal") 
      opt->SetMscStepLimitation(fMinimal);

    else if(newValue == "UseDistanceToBoundary") 
      opt->SetMscStepLimitation(fUseDistanceToBoundary);

    else if(newValue == "UseSafety")
      opt->SetMscStepLimitation(fUseSafety);

    else {
      G4cout << "### G4EnergyLossMessenger WARNING: StepLimit type <" 
	     << newValue << "> unknown!" << G4endl;
      return;
    }
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == MinEnCmd) {
    opt->SetMinEnergy(MinEnCmd->GetNewDoubleValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == MaxEnCmd) { 
    opt->SetMaxEnergy(MaxEnCmd->GetNewDoubleValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }  
  if (command == IntegCmd) {
    opt->SetIntegral(IntegCmd->GetNewBoolValue(newValue));
    return;
  }
  if (command == rangeCmd) {
    opt->SetBuildCSDARange(rangeCmd->GetNewBoolValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }  
  if (command == lpmCmd) {
    opt->SetLPMFlag(lpmCmd->GetNewBoolValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == splCmd) {
    opt->SetSplineFlag(splCmd->GetNewBoolValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == aplCmd) {
    opt->SetApplyCuts(aplCmd->GetNewBoolValue(newValue));
    return;
  }
  if (command == latCmd) {
    opt->SetMscLateralDisplacement(latCmd->GetNewBoolValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == verCmd) {
    opt->SetVerbose(verCmd->GetNewIntValue(newValue));
    return;
  }
  if (command == ver1Cmd) {
    opt->SetVerbose(ver1Cmd->GetNewIntValue(newValue));
    return;
  }
  if (command == lllCmd) {
    opt->SetLinearLossLimit(lllCmd->GetNewDoubleValue(newValue));
    return;
  }
  if (command == labCmd) {
    opt->SetLambdaFactor(labCmd->GetNewDoubleValue(newValue));
    return;
  }
  if (command == skinCmd) {
    opt->SetSkin(skinCmd->GetNewDoubleValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }  
  if (command == dedxCmd) { 
    opt->SetDEDXBinning(dedxCmd->GetNewIntValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == lamCmd) { 
    opt->SetLambdaBinning(lamCmd->GetNewIntValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == frCmd) {
    opt->SetMscRangeFactor(frCmd->GetNewDoubleValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == fgCmd) {
    opt->SetMscGeomFactor(fgCmd->GetNewDoubleValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == mscfCmd) {
    opt->SetFactorForAngleLimit(mscfCmd->GetNewDoubleValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    return;
  }
  if (command == angCmd) { 
    opt->SetPolarAngleLimit(angCmd->GetNewDoubleValue(newValue));
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
