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
// File name:     G4EmExtraParametersMessenger
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 07-05-2019 
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmExtraParametersMessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImanager.hh"
#include "G4EmExtraParameters.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmExtraParametersMessenger::G4EmExtraParametersMessenger(G4EmExtraParameters* ptr) 
  : theParameters(ptr)
{
  paiCmd = new G4UIcommand("/process/em/AddPAIRegion",this);
  paiCmd->SetGuidance("Activate PAI in the G4Region.");
  paiCmd->SetGuidance("  partName  : particle name (default - all)");
  paiCmd->SetGuidance("  regName   : G4Region name");
  paiCmd->SetGuidance("  paiType   : PAI, PAIphoton");
  paiCmd->AvailableForStates(G4State_PreInit);

  G4UIparameter* part = new G4UIparameter("partName",'s',false);
  paiCmd->SetParameter(part);

  G4UIparameter* pregName = new G4UIparameter("regName",'s',false);
  paiCmd->SetParameter(pregName);

  G4UIparameter* ptype = new G4UIparameter("type",'s',false);
  paiCmd->SetParameter(ptype);
  ptype->SetParameterCandidates("pai PAI PAIphoton");

  mscoCmd = new G4UIcommand("/process/em/AddEmRegion",this);
  mscoCmd->SetGuidance("Add optional EM configuration for a G4Region.");
  mscoCmd->SetGuidance("  regName  : G4Region name");
  mscoCmd->SetGuidance("  emType   : G4EmStandard, G4EmStandard_opt1, ...");
  mscoCmd->AvailableForStates(G4State_PreInit);

  G4UIparameter* mregName = new G4UIparameter("regName",'s',false);
  mscoCmd->SetParameter(mregName);

  G4UIparameter* mtype = new G4UIparameter("mscType",'s',false);
  mscoCmd->SetParameter(mtype);
  mtype->SetParameterCandidates("G4EmStandard G4EmStandard_opt1 G4EmStandard_opt2 G4EmStandard_opt3 G4EmStandard_opt4 G4EmStandardGS G4EmStandardSS G4EmLivermore G4EmPenelope G4RadioactiveDecay");

  SubSecCmd = new G4UIcommand("/process/eLoss/subsec",this);
  SubSecCmd->SetGuidance("Switch true/false the subcutoff generation per region.");
  SubSecCmd->SetGuidance("  subSec   : true/false");
  SubSecCmd->SetGuidance("  Region   : region name");
  SubSecCmd->AvailableForStates(G4State_PreInit);

  G4UIparameter* subSec = new G4UIparameter("subSec",'s',false);
  SubSecCmd->SetParameter(subSec);

  G4UIparameter* subSecReg = new G4UIparameter("Region",'s',false);
  SubSecCmd->SetParameter(subSecReg);

  StepFuncCmd = new G4UIcommand("/process/eLoss/StepFunction",this);
  StepFuncCmd->SetGuidance("Set the energy loss step limitation parameters for e+-.");
  StepFuncCmd->SetGuidance("  dRoverR   : max Range variation per step");
  StepFuncCmd->SetGuidance("  finalRange: range for final step");
  StepFuncCmd->SetGuidance("  unit      : unit of finalRange");
  StepFuncCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* dRoverRPrm = new G4UIparameter("dRoverR",'d',false);
  dRoverRPrm->SetParameterRange("dRoverR>0. && dRoverR<=1.");
  StepFuncCmd->SetParameter(dRoverRPrm);

  G4UIparameter* finalRangePrm = new G4UIparameter("finalRange",'d',false);
  finalRangePrm->SetParameterRange("finalRange>0.");
  StepFuncCmd->SetParameter(finalRangePrm);

  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',true);
  unitPrm->SetDefaultUnit("mm");
  StepFuncCmd->SetParameter(unitPrm);

  StepFuncCmd1 = new G4UIcommand("/process/eLoss/StepFunctionMuHad",this);
  StepFuncCmd1->SetGuidance("Set the energy loss step limitation parameters for muon/hadron.");
  StepFuncCmd1->SetGuidance("  dRoverR   : max Range variation per step");
  StepFuncCmd1->SetGuidance("  finalRange: range for final step");
  StepFuncCmd1->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* dRoverRPrm1 = new G4UIparameter("dRoverRMuHad",'d',false);
  dRoverRPrm1->SetParameterRange("dRoverRMuHad>0. && dRoverRMuHad<=1.");
  StepFuncCmd1->SetParameter(dRoverRPrm1);

  G4UIparameter* finalRangePrm1 = new G4UIparameter("finalRangeMuHad",'d',false);
  finalRangePrm1->SetParameterRange("finalRangeMuHad>0.");
  StepFuncCmd1->SetParameter(finalRangePrm1);

  G4UIparameter* unitPrm1 = new G4UIparameter("unit",'s',true);
  unitPrm1->SetDefaultValue("mm");
  StepFuncCmd1->SetParameter(unitPrm1);

  StepFuncCmd2 = new G4UIcommand("/process/eLoss/StepFunctionLightIons",this);
  StepFuncCmd2->SetGuidance("Set the energy loss step limitation parameters for light ions.");
  StepFuncCmd2->SetGuidance("  dRoverR   : max Range variation per step");
  StepFuncCmd2->SetGuidance("  finalRange: range for final step");
  StepFuncCmd2->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* dRoverRPrm2 = new G4UIparameter("dRoverRLIons",'d',false);
  dRoverRPrm2->SetParameterRange("dRoverRLIons>0. && dRoverRLIons<=1.");
  StepFuncCmd2->SetParameter(dRoverRPrm2);

  G4UIparameter* finalRangePrm2 = new G4UIparameter("finalRangeLIons",'d',false);
  finalRangePrm2->SetParameterRange("finalRangeLIons>0.");
  StepFuncCmd2->SetParameter(finalRangePrm2);

  G4UIparameter* unitPrm2 = new G4UIparameter("unit",'s',true);
  unitPrm2->SetDefaultValue("mm");
  StepFuncCmd2->SetParameter(unitPrm2);

  StepFuncCmd3 = new G4UIcommand("/process/eLoss/StepFunctionIons",this);
  StepFuncCmd3->SetGuidance("Set the energy loss step limitation parameters for ions.");
  StepFuncCmd3->SetGuidance("  dRoverR   : max Range variation per step");
  StepFuncCmd3->SetGuidance("  finalRange: range for final step");
  StepFuncCmd3->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* dRoverRPrm3 = new G4UIparameter("dRoverRMuHad",'d',false);
  dRoverRPrm3->SetParameterRange("dRoverRIons>0. && dRoverRIons<=1.");
  StepFuncCmd3->SetParameter(dRoverRPrm3);

  G4UIparameter* finalRangePrm3 = new G4UIparameter("finalRangeIons",'d',false);
  finalRangePrm3->SetParameterRange("finalRangeIons>0.");
  StepFuncCmd3->SetParameter(finalRangePrm3);

  G4UIparameter* unitPrm3 = new G4UIparameter("unit",'s',true);
  unitPrm3->SetDefaultValue("mm");
  StepFuncCmd3->SetParameter(unitPrm3);

  bfCmd = new G4UIcommand("/process/em/setBiasingFactor",this);
  bfCmd->SetGuidance("Set factor for the process cross section.");
  bfCmd->SetGuidance("  procName   : process name");
  bfCmd->SetGuidance("  procFact   : factor");
  bfCmd->SetGuidance("  flagFact   : flag to change weight");
  bfCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* procName = new G4UIparameter("procName",'s',false);
  bfCmd->SetParameter(procName);

  G4UIparameter* procFact = new G4UIparameter("procFact",'d',false);
  bfCmd->SetParameter(procFact);

  G4UIparameter* flagFact = new G4UIparameter("flagFact",'s',false);
  bfCmd->SetParameter(flagFact);

  fiCmd = new G4UIcommand("/process/em/setForcedInteraction",this);
  fiCmd->SetGuidance("Set factor for the process cross section.");
  fiCmd->SetGuidance("  procNam    : process name");
  fiCmd->SetGuidance("  regNam     : region name");
  fiCmd->SetGuidance("  tlength    : fixed target length");
  fiCmd->SetGuidance("  unitT      : length unit");
  fiCmd->SetGuidance("  tflag      : flag to change weight");
  fiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* procNam = new G4UIparameter("procNam",'s',false);
  fiCmd->SetParameter(procNam);

  G4UIparameter* regNam = new G4UIparameter("regNam",'s',false);
  fiCmd->SetParameter(regNam);

  G4UIparameter* tlength = new G4UIparameter("tlength",'d',false);
  tlength->SetParameterRange("tlength>0");
  fiCmd->SetParameter(tlength);

  G4UIparameter* unitT = new G4UIparameter("unitT",'s',true);
  unitT->SetDefaultUnit("mm");
  fiCmd->SetParameter(unitT);

  G4UIparameter* flagT = new G4UIparameter("tflag",'b',true);
  flagT->SetDefaultValue(true);
  fiCmd->SetParameter(flagT);

  bsCmd = new G4UIcommand("/process/em/setSecBiasing",this);
  bsCmd->SetGuidance("Set bremsstrahlung or delta-e- splitting/Russian roulette per region.");
  bsCmd->SetGuidance("  bProcNam : process name");
  bsCmd->SetGuidance("  bRegNam  : region name");
  bsCmd->SetGuidance("  bFactor  : number of split gamma or probability of Russian roulette");
  bsCmd->SetGuidance("  bEnergy  : max energy of a secondary for this biasing method");
  bsCmd->SetGuidance("  bUnit    : energy unit");
  bsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  G4UIparameter* bProcNam = new G4UIparameter("bProcNam",'s',false);
  bsCmd->SetParameter(bProcNam);

  G4UIparameter* bRegNam = new G4UIparameter("bRegNam",'s',false);
  bsCmd->SetParameter(bRegNam);

  G4UIparameter* bFactor = new G4UIparameter("bFactor",'d',false);
  bsCmd->SetParameter(bFactor);

  G4UIparameter* bEnergy = new G4UIparameter("bEnergy",'d',false);
  bsCmd->SetParameter(bEnergy);

  G4UIparameter* bUnit = new G4UIparameter("bUnit",'s',true);
  bUnit->SetDefaultUnit("MeV");
  bsCmd->SetParameter(bUnit);

  dirSplitCmd = new G4UIcmdWithABool("/process/em/setDirectionalSplitting",this);
  dirSplitCmd->SetGuidance("Enable directional brem splitting");
  dirSplitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  qeCmd = new G4UIcmdWithABool("/process/em/QuantumEntanglement",this);
  qeCmd->SetGuidance("Enable quantum entanglement");
  qeCmd->AvailableForStates(G4State_PreInit);

  dirSplitTargetCmd = new G4UIcmdWith3VectorAndUnit("/process/em/setDirectionalSplittingTarget",this);
  dirSplitTargetCmd->SetGuidance("Position of arget for directional splitting");
  dirSplitTargetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  dirSplitRadiusCmd = new G4UIcmdWithADoubleAndUnit("/process/em/setDirectionalSplittingRadius",this);
  dirSplitRadiusCmd->SetGuidance("Radius of target for directional splitting");
  dirSplitRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmExtraParametersMessenger::~G4EmExtraParametersMessenger()
{
  delete paiCmd;
  delete mscoCmd;
  delete SubSecCmd;
  delete bfCmd;
  delete fiCmd;
  delete bsCmd;
  delete qeCmd;
  delete StepFuncCmd;
  delete StepFuncCmd1;
  delete StepFuncCmd2;
  delete StepFuncCmd3;
  delete dirSplitCmd;
  delete dirSplitTargetCmd;
  delete dirSplitRadiusCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmExtraParametersMessenger::SetNewValue(G4UIcommand* command, 
                                               G4String newValue)
{
  G4bool physicsModified = false;

  if (command == paiCmd) {
    G4String s1(""),s2(""),s3("");
    std::istringstream is(newValue);
    is >> s1 >> s2 >> s3;
    theParameters->AddPAIModel(s1, s2, s3);
  } else if (command == mscoCmd) {
    G4String s1(""),s2("");
    std::istringstream is(newValue);
    is >> s1 >> s2;
    theParameters->AddPhysics(s1, s2);
  } else if (command == StepFuncCmd || command == StepFuncCmd1 || command == StepFuncCmd2 || command == StepFuncCmd3) {
    G4double v1,v2;
    G4String unt;
    std::istringstream is(newValue);
    is >> v1 >> v2 >> unt;
    v2 *= G4UIcommand::ValueOf(unt);
    if(command == StepFuncCmd) {
      theParameters->SetStepFunction(v1,v2);
    } else if(command == StepFuncCmd1) {
      theParameters->SetStepFunctionMuHad(v1,v2);
    } else if(command == StepFuncCmd2) {
      theParameters->SetStepFunctionLightIons(v1,v2);
    } else {
      theParameters->SetStepFunctionIons(v1,v2);
    }
    physicsModified = true;
  } else if (command == SubSecCmd) {
    G4String s1, s2;
    std::istringstream is(newValue);
    is >> s1 >> s2;
    G4bool yes = false;
    if(s1 == "true") { yes = true; }
    theParameters->SetSubCutoff(yes,s2);
  } else if (command == bfCmd) {
    G4double v1(1.0);
    G4String s0(""),s1("");
    std::istringstream is(newValue);
    is >> s0 >> v1 >> s1;
    G4bool yes = false;
    if(s1 == "true") { yes = true; }
    theParameters->SetProcessBiasingFactor(s0,v1,yes);
    physicsModified = true;
  } else if (command == fiCmd) {
    G4double v1(0.0);
    G4String s1(""),s2(""),s3(""),unt("mm");
    std::istringstream is(newValue);
    is >> s1 >> s2 >> v1 >> unt >> s3;
    G4bool yes = false;
    if(s3 == "true") { yes = true; }
    v1 *= G4UIcommand::ValueOf(unt);
    theParameters->ActivateForcedInteraction(s1,s2,v1,yes);
    physicsModified = true;
  } else if (command == bsCmd) {
    G4double fb(1.0),en(1.e+30);
    G4String s1(""),s2(""),unt("MeV");
    std::istringstream is(newValue);
    is >> s1 >> s2 >> fb >> en >> unt;
    en *= G4UIcommand::ValueOf(unt);    
    theParameters->ActivateSecondaryBiasing(s1,s2,fb,en);
    physicsModified = true;
  } else if (command == qeCmd) {
    theParameters->SetQuantumEntanglement(qeCmd->GetNewBoolValue(newValue));
  } else if (command == dirSplitCmd) {
    theParameters->SetDirectionalSplitting(
      dirSplitCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == dirSplitTargetCmd) {
    G4ThreeVector t = dirSplitTargetCmd->GetNew3VectorValue(newValue);
    theParameters->SetDirectionalSplittingTarget(t);
    physicsModified = true;
  } else if (command == dirSplitRadiusCmd) {
    G4double r = dirSplitRadiusCmd->GetNewDoubleValue(newValue);
    theParameters->SetDirectionalSplittingRadius(r);
    physicsModified = true;
  }
  
  if(physicsModified) {
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
