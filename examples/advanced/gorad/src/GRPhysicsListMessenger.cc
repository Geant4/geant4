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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRPhysicsListMessenger.cc
//   A messenger class that handles Gorad physics list options.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRPhysicsListMessenger.hh"

#include "GRPhysicsList.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

GRPhysicsListMessenger::GRPhysicsListMessenger(GRPhysicsList* pl)
: pPL(pl)
{
  G4UIparameter* param = nullptr;

  physDir = new G4UIdirectory("/gorad/physics/");
  physDir->SetGuidance("GORAD physics selection");

  selectEMCmd = new G4UIcmdWithAString("/gorad/physics/EM",this);
  selectEMCmd->AvailableForStates(G4State_PreInit);
  selectEMCmd->SetToBeBroadcasted(false);
  selectEMCmd->SetParameterName("EM_option",true);
  selectEMCmd->SetCandidates("Op_0 Op_1 Op_3 Op_4 LIV LIV_Pol");
  selectEMCmd->SetDefaultValue("Op_0");
  selectEMCmd->SetGuidance("Select EM Physics option");
  selectEMCmd->SetGuidance(" Op_0 (default) : Suitable to medium and high energy applications");
  selectEMCmd->SetGuidance(" Op_1 : Faster than Op_0 because of less accurate MSC step limitation");
  selectEMCmd->SetGuidance(" Op_3 : Suitable for medical applications - more accurate MSC for all particles");
  selectEMCmd->SetGuidance(" Op_4 : Most accurate (GS MSC model with Mott correction and error-free stepping for e+/-");
  selectEMCmd->SetGuidance(" LIV  : Livermore models for e-/gamma below 1 GeV, otherwise Op_0");
  selectEMCmd->SetGuidance(" LIV_Pol : Polarized extension of Livermore models (t.b.a.)");

  selectHadCmd = new G4UIcmdWithAString("/gorad/physics/Hadronic",this);
  selectHadCmd->AvailableForStates(G4State_PreInit);
  selectHadCmd->SetToBeBroadcasted(false);
  selectHadCmd->SetParameterName("Had_option",true);
  selectHadCmd->SetCandidates("FTFP_BERT QGSP_BIC Shielding");
  selectHadCmd->SetDefaultValue("FTFP_BERT");
  selectHadCmd->SetGuidance("Select Hadronic Physics option");
  selectHadCmd->SetGuidance(" FTFP_BERT (default) : Fritiof string + Bertini cascade + Precompound de-excitation");
  selectHadCmd->SetGuidance("                       suitable to most of midium and high energy applications");
  selectHadCmd->SetGuidance(" QGSP_BIC : Quark-Gluon-String + Fritiof string + Binary cascade + Precompound de-excitation");
  selectHadCmd->SetGuidance("            suitable for lower energy applications such as medical");
  selectHadCmd->SetGuidance(" Shielding : Similar to FTFP+BERT with better ion-ion interactions.");
  selectHadCmd->SetGuidance("             High-Precision neutron and Radioactive Decay models are included by default.");

  addHPCmd = new G4UIcmdWithoutParameter("/gorad/physics/addHP",this);
  addHPCmd->AvailableForStates(G4State_PreInit);
  addHPCmd->SetToBeBroadcasted(false);
  addHPCmd->SetGuidance("Add High-Precision neutron model.");
  addHPCmd->SetGuidance(" Note: Shielding option has already had HP. This command does not make effect to Shielding option.");

  addRDMCmd = new G4UIcmdWithoutParameter("/gorad/physics/addRDM",this);
  addRDMCmd->AvailableForStates(G4State_PreInit);
  addRDMCmd->SetToBeBroadcasted(false);
  addRDMCmd->SetGuidance("Add Radioactive Decay model.");
  addRDMCmd->SetGuidance(" Note: Shielding option has already had RDM. This command does not make effect to Shielding option.");

  addRMCCmd = new G4UIcmdWithoutParameter("/gorad/physics/addRMC",this);
  addRMCCmd->AvailableForStates(G4State_PreInit);
  addRMCCmd->SetToBeBroadcasted(false);
  addRMCCmd->SetGuidance("Add Reverse Monte Carlo.");

  addOpticalCmd = new G4UIcmdWithoutParameter("/gorad/physics/addOptical",this);
  addOpticalCmd->AvailableForStates(G4State_PreInit);
  addOpticalCmd->SetToBeBroadcasted(false);
  addOpticalCmd->SetGuidance("Add Optical physics");

  addStepLimitCmd = new G4UIcmdWithAString("/gorad/physics/addStepLimit",this);
  addStepLimitCmd->AvailableForStates(G4State_PreInit);
  addStepLimitCmd->SetToBeBroadcasted(false);
  addStepLimitCmd->SetGuidance("Add step-limiter process to artificially limit step length.");
  addStepLimitCmd->SetGuidance("Specify particle types to be applied.");
  addStepLimitCmd->SetGuidance("  charged (default) : applied only to the charged particles");
  addStepLimitCmd->SetGuidance("  neutral : applied only to the neutral particles");
  addStepLimitCmd->SetGuidance("  all : applied to all particle types");
  addStepLimitCmd->SetGuidance("  e+/- : applied only to e+/e-");
  addStepLimitCmd->SetGuidance(" Note: In addition to this command, you need to specify the limitation value by");
  addStepLimitCmd->SetGuidance("       /gorad/physics/limit/stepLimit or /gorad/physics/limit/localStepLimt command.");
  addStepLimitCmd->SetParameterName("particle",true);
  addStepLimitCmd->SetDefaultValue("charged");
  addStepLimitCmd->SetCandidates("charged neutral all e+/-");

  physLimitDir = new G4UIdirectory("/gorad/physics/limit/");
  physLimitDir->SetGuidance("Specify step limitation");

  setStepLimitCmd = new G4UIcmdWithADoubleAndUnit("/gorad/physics/limit/stepLimit",this);
  setStepLimitCmd->AvailableForStates(G4State_Idle);
  setStepLimitCmd->SetToBeBroadcasted(false);
  setStepLimitCmd->SetParameterName("length",false);
  setStepLimitCmd->SetDefaultUnit("mm");
  setStepLimitCmd->SetGuidance("Define the limitation of the step length");
  setStepLimitCmd->SetGuidance("This limitation is applied to the entire geometry except regions that has its dedicated limit.");

  setRegionStepLimitCmd = new G4UIcommand("/gorad/physics/limit/regionStepLimit",this);
  setRegionStepLimitCmd->AvailableForStates(G4State_Idle);
  setRegionStepLimitCmd->SetToBeBroadcasted(false);
  setRegionStepLimitCmd->SetGuidance("Define the limitation of the step length for the specified region");
  setRegionStepLimitCmd->SetGuidance("   [usage] /gorad/physics/limit/regionStepLimit region length [unit]");
  setRegionStepLimitCmd->SetGuidance("      region (string) : region name");
  setRegionStepLimitCmd->SetGuidance(" Note: Region has to be defined in advance to this command.");
  setRegionStepLimitCmd->SetGuidance("       If new region is necessary, use /gorad/geometry/createRegion to create it.");
  param = new G4UIparameter("region",'s',false);
  setRegionStepLimitCmd->SetParameter(param);
  param = new G4UIparameter("length",'d',false);
  setRegionStepLimitCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("mm");
  setRegionStepLimitCmd->SetParameter(param);

  physCutDir = new G4UIdirectory("/gorad/physics/cuts/");
  physCutDir->SetGuidance("Specify production thresholds (a.k.a. cuts)");

  setCutCmd = new G4UIcmdWithADoubleAndUnit("/gorad/physics/cuts/setCuts",this);
  setCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  setCutCmd->SetToBeBroadcasted(false);
  setCutCmd->SetParameterName("length",false);
  setCutCmd->SetDefaultUnit("mm");
  setCutCmd->SetGuidance("Specify production thresholds (a.k.a. cuts) that is applied to the entire geometry");
  setCutCmd->SetGuidance("This threshold is applied to all of e-, e+, gamma and proton.");
  setCutCmd->SetGuidance("Threshold of each particle can be overwitted by /gorad/physics/cuts/setParticleCut command");

  setCutParticleCmd = new G4UIcommand("/gorad/physics/cuts/setParticleCut",this);
  setCutParticleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  setCutParticleCmd->SetToBeBroadcasted(false);
  setCutParticleCmd->SetGuidance("Specify production threshold (a.k.a. cut) for the specified particle that is applied to the entire geometry");
  setCutParticleCmd->SetGuidance("  [usage] /gorad/physics/setParticleCut particle cut unit");
  param = new G4UIparameter("particle",'s',false);
  param->SetParameterCandidates("e- e+ gamma proton");
  setCutParticleCmd->SetParameter(param);
  param = new G4UIparameter("cut",'d',false);
  setCutParticleCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("mm");
  setCutParticleCmd->SetParameter(param);

  setCutRegionCmd = new G4UIcommand("/gorad/physics/cuts/setRegionCut",this);
  setCutRegionCmd->AvailableForStates(G4State_Idle);
  setCutRegionCmd->SetToBeBroadcasted(false);
  setCutRegionCmd->SetGuidance("Specify production threshold (a.k.a. cut) that is applied to the specified region");
  setCutRegionCmd->SetGuidance("  [usage] /gorad/physics/setRegionCut region cut unit");
  setCutRegionCmd->SetGuidance("This threshold is applied to all of e-, e+, gamma and proton.");
  setCutRegionCmd->SetGuidance("Threshold of each particle can be overwitted by /gorad/physics/cuts/setRegionParticleCut command");
  setCutRegionCmd->SetGuidance(" Note: Region has to be defined in advance to this command.");
  setCutRegionCmd->SetGuidance("       If new region is necessary, use /gorad/geometry/createRegion to create it.");
  param = new G4UIparameter("region",'s',false);
  setCutRegionCmd->SetParameter(param);
  param = new G4UIparameter("cut",'d',false);
  setCutRegionCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("mm");
  setCutRegionCmd->SetParameter(param);

  setCutRegionParticleCmd = new G4UIcommand("/gorad/physics/cuts/setRegionParticleCut",this);
  setCutRegionParticleCmd->AvailableForStates(G4State_Idle);
  setCutRegionParticleCmd->SetToBeBroadcasted(false);
  setCutRegionParticleCmd->SetGuidance("Specify production threshold (a.k.a. cut) that is applied to the specified region");
  setCutRegionParticleCmd->SetGuidance("  [usage] /gorad/physics/setRegionParticleCut region particle cut unit");
  setCutRegionParticleCmd->SetGuidance(" Note: Region has to be defined in advance to this command.");
  setCutRegionParticleCmd->SetGuidance("       If new region is necessary, use /gorad/geometry/createRegion to create it.");
  param = new G4UIparameter("region",'s',false);
  setCutRegionParticleCmd->SetParameter(param);
  param = new G4UIparameter("particle",'s',false);
  param->SetParameterCandidates("e- e+ gamma proton");
  setCutRegionParticleCmd->SetParameter(param);
  param = new G4UIparameter("cut",'d',false);
  setCutRegionParticleCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("mm");
  setCutRegionParticleCmd->SetParameter(param);

}

GRPhysicsListMessenger::~GRPhysicsListMessenger()
{
  delete selectEMCmd;
  delete selectHadCmd;
  delete addHPCmd;
  delete addRDMCmd;
  delete addRMCCmd;
  delete addOpticalCmd;
  delete addStepLimitCmd;
  delete setStepLimitCmd;
  delete setRegionStepLimitCmd;
  delete setCutCmd;
  delete setCutParticleCmd;
  delete setCutRegionCmd;
  delete setCutRegionParticleCmd;

  delete physLimitDir;
  delete physCutDir;
  delete physDir;
}

#include "G4Tokenizer.hh"

void GRPhysicsListMessenger::SetNewValue(G4UIcommand* cmd, G4String val)
{
  if(cmd==selectEMCmd)
  { pPL->SetEM(val); }
  else if(cmd==selectHadCmd)
  { pPL->SetHad(val); }
  else if(cmd==addHPCmd)
  { pPL->AddHP(); }
  else if(cmd==addRDMCmd)
  { pPL->AddRDM(); }
  else if(cmd==addRMCCmd)
  { pPL->AddRMC(); }
  else if(cmd==addOpticalCmd)
  { G4cout<<"Not yet implemented."<<G4endl; }
  else if(cmd==addStepLimitCmd)
  {
    G4int opt = 0;
    if(val=="neutral") opt = 1; 
    else if(val=="all") opt = 2; 
    else if(val=="e+/-") opt = 3; 
    pPL->AddStepLimit(opt);
  }
  else if(cmd==setStepLimitCmd)
  { pPL->SetGlobalStepLimit(setStepLimitCmd->GetNewDoubleValue(val)); }
  else if(cmd==setRegionStepLimitCmd)
  {
    G4Tokenizer next(val);
    G4String reg = next();
    G4String newVal = next();
    newVal += " ";
    newVal += next();
    auto regPtr = pPL->SetLocalStepLimit(reg,setRegionStepLimitCmd->ConvertToDimensionedDouble(newVal));
    if(!regPtr)
    {
      G4ExceptionDescription ed;
      ed << "Region <" << reg << "> is not defined. Region has to be defined in advance to this command."
         << "\nIf new region is necessary, use /gorad/geometry/createRegion to create it.";
      setRegionStepLimitCmd->CommandFailed(ed);
    }
  }
  else if(cmd==setCutCmd)
  { pPL->SetGlobalCuts(setCutCmd->GetNewDoubleValue(val)); }
  else if(cmd==setCutParticleCmd)
  {
    G4Tokenizer next(val);
    G4String pat = next();
    G4String newVal = next();
    newVal += " ";
    newVal += next();
    G4int i = 0;
    if(pat=="e-") i = 0; 
    else if(pat=="e+") i = 1; 
    else if(pat=="gamma") i = 2; 
    else if(pat=="proton") i = 3; 
    pPL->SetGlobalCut(i,setCutParticleCmd->ConvertToDimensionedDouble(newVal));
  }
  else if(cmd==setCutRegionCmd)
  {
    G4Tokenizer next(val);
    G4String reg = next();
    G4String newVal = next();
    newVal += " ";
    newVal += next();
    auto regPtr = pPL->SetLocalCuts(reg,setCutRegionCmd->ConvertToDimensionedDouble(newVal));
    if(!regPtr)
    {
      G4ExceptionDescription ed;
      ed << "Region <" << reg << "> is not defined. Region has to be defined in advance to this command."
         << "\nIf new region is necessary, use /gorad/geometry/createRegion to create it.";
      setRegionStepLimitCmd->CommandFailed(ed);
    }
  }
  else if(cmd==setCutRegionParticleCmd)
  {
    G4Tokenizer next(val);
    G4String reg = next();
    G4String pat = next();
    G4int i = 0;
    if(pat=="e-") i = 0; 
    else if(pat=="e+") i = 1; 
    else if(pat=="gamma") i = 2; 
    else if(pat=="proton") i = 3; 
    G4String newVal = next();
    newVal += " ";
    newVal += next();
    auto regPtr = pPL->SetLocalCut(reg,i,setCutRegionParticleCmd->ConvertToDimensionedDouble(newVal));
    if(!regPtr)
    {
      G4ExceptionDescription ed;
      ed << "Region <" << reg << "> is not defined. Region has to be defined in advance to this command."
         << "\nIf new region is necessary, use /gorad/geometry/createRegion to create it.";
      setRegionStepLimitCmd->CommandFailed(ed);
    }
  }
  
}

G4String GRPhysicsListMessenger::GetCurrentValue(G4UIcommand* cmd)
{
  G4String val("");

  if(cmd==selectEMCmd)
  { val = pPL->GetEM(); }
  else if(cmd==selectHadCmd)
  { val = pPL->GetHad(); }
  else if(cmd==addHPCmd)
  { val = cmd->ConvertToString(pPL->IfHP()); }
  else if(cmd==addRDMCmd)
  { val = cmd->ConvertToString(pPL->IfRDM()); }
  else if(cmd==addRMCCmd)
  { val = cmd->ConvertToString(pPL->IfRMC()); }
  else if(cmd==addOpticalCmd)
  { G4cout<<"Not yet implemented."<<G4endl; }
  else if(cmd==addStepLimitCmd)
  {
    auto opt = pPL->IfStepLimit();
    switch(opt)
    {
      case 0: val =  "charged"; break;
      case 1: val =  "neutral"; break;
      case 2: val =  "all"; break;
      case 3: val =  "e+/-"; break;
      default : val = "undefined"; break;
    }
  }
  return val;
}


