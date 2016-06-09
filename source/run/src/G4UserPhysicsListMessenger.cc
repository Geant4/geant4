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
// $Id: G4UserPhysicsListMessenger.cc,v 1.24 2004/11/01 16:47:56 kurasige Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// 
//---------------------------------------------------------------
//
//  G4UserPhysicsListMessenger.cc
// ------------------------------------------------------------
//	History
//        first version                   09 Jan. 1998 by H.Kurashige 
//        add buildPhysicsTable command   13 Apr. 1999 by H.Kurashige
//        add setStoredInAscii command    12 Mar. 2001 by H.Kurashige
// ------------------------------------------------------------

#include "G4UserPhysicsListMessenger.hh"
#include "G4VUserPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"
#include "G4Tokenizer.hh"           


G4UserPhysicsListMessenger::G4UserPhysicsListMessenger(G4VUserPhysicsList* pParticleList):thePhysicsList(pParticleList)
{
  G4UIparameter* param = 0;
  // /run/particle    directory
  theDirectory = new G4UIdirectory("/run/particle/");
  theDirectory->SetGuidance("Commands for G4VUserPhysicsList.");

  // /run/particle/Verbose command
  verboseCmd = new G4UIcmdWithAnInteger("/run/particle/verbose",this);
  verboseCmd->SetGuidance("Set the Verbose level of G4VUserPhysicsList.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display warning messages");
  verboseCmd->SetGuidance(" 2 : Display more");
  verboseCmd->SetParameterName("level",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("level >=0 && level <=3");

  // /run/particle/setCut command
  setPCutCmd = new G4UIcmdWithADoubleAndUnit("/run/particle/setCut",this);
  setPCutCmd->SetGuidance("Set default cut value ");
  setPCutCmd->SetGuidance("This command is equivallent to /run/setCut command.");
  setPCutCmd->SetGuidance("This command is kept for backward compatibility.");
  setPCutCmd->SetParameterName("cut",false);
  setPCutCmd->SetDefaultValue(1.0);
  setPCutCmd->SetRange("cut >0.0");
  setPCutCmd->SetDefaultUnit("mm");
  setPCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // /run/setCut command
  setCutCmd = new G4UIcmdWithADoubleAndUnit("/run/setCut",this);
  setCutCmd->SetGuidance("Set default cut value ");
  setCutCmd->SetParameterName("cut",false);
  setCutCmd->SetDefaultValue(1.0);
  setCutCmd->SetRange("cut >0.0");
  setCutCmd->SetDefaultUnit("mm");
  setCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // /run/setCutForRegion command
  setCutRCmd = new G4UIcommand("/run/setCutForRegion",this);
  setCutRCmd->SetGuidance("Set cut value for a region");
  param = new G4UIparameter("Region",'s',false);
  setCutRCmd->SetParameter(param);
  param = new G4UIparameter("cut",'d',false);
  param->SetParameterRange("cut >0.0");
  setCutRCmd->SetParameter(param);
  param = new G4UIparameter("Unit",'s',true);
  param->SetDefaultValue("mm");
  param->SetParameterCandidates(setCutRCmd->UnitsList(setCutRCmd->CategoryOf("mm")));
  setCutRCmd->SetParameter(param);
  setCutRCmd->AvailableForStates(G4State_Idle);

  // /run/particle/DumpList command
  dumpListCmd = new G4UIcmdWithoutParameter("/run/particle/dumpList",this);
  dumpListCmd->SetGuidance("Dump List of particles in G4VUserPhysicsList. ");

  // /run/particle/addProcManager command
  addProcManCmd = new G4UIcmdWithAString("/run/particle/addProcManager", this);
  addProcManCmd->SetGuidance("add process manager to specified particle type");
  addProcManCmd->SetParameterName("particleType", true);
  addProcManCmd->SetDefaultValue("");
  addProcManCmd->AvailableForStates(G4State_Init,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  // /run/particle/buildPhysicsTable command
  buildPTCmd = new G4UIcmdWithAString("/run/particle/buildPhysicsTable", this);
  buildPTCmd->SetGuidance("build physics table of specified particle type");
  buildPTCmd->SetParameterName("particleType", true);
  buildPTCmd->SetDefaultValue("");
  buildPTCmd->AvailableForStates(G4State_Init,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  // /run/particle/storePhysicsTable command
  storeCmd = new G4UIcmdWithAString("/run/particle/storePhysicsTable",this);
  storeCmd->SetGuidance("Store Physics Table");
  storeCmd->SetGuidance("  Enter directory name");
  storeCmd->SetParameterName("dirName",true);
  storeCmd->SetDefaultValue("");
  storeCmd->AvailableForStates(G4State_Idle);

  //  /run/particle/retrievePhysicsTable command
  retrieveCmd = new G4UIcmdWithAString("/run/particle/retrievePhysicsTable",this);
  retrieveCmd->SetGuidance("Retrieve Physics Table");
  retrieveCmd->SetGuidance("  Enter directory name or OFF to switch off");
  retrieveCmd->SetParameterName("dirName",true);
  retrieveCmd->SetDefaultValue("");
  retrieveCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //  /run/particle/setStoredInAscii command
  asciiCmd = new G4UIcmdWithAnInteger("/run/particle/setStoredInAscii",this);
  asciiCmd->SetGuidance("Switch on/off ascii mode in store/retreive Physics Table");
  asciiCmd->SetGuidance("  Enter 0(binary) or 1(ascii)");
  asciiCmd->SetParameterName("ascii",true);
  asciiCmd->SetDefaultValue(0);
  asciiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  asciiCmd->SetRange("ascii ==0 || ascii ==1");

  //Commnad    /run/particle/applyCuts command
  applyCutsCmd = new G4UIcommand("/run/particle/applyCuts",this);
  applyCutsCmd->SetGuidance("Set applyCuts flag for a particle.");
  applyCutsCmd->SetGuidance(" Some EM processes which do not have infrared divergence");
  applyCutsCmd->SetGuidance("may generate gamma, e- and/or e+ with kinetic energies");
  applyCutsCmd->SetGuidance("below the production threshold. By setting this flag,");
  applyCutsCmd->SetGuidance("such secondaries below threshold are eliminated and");
  applyCutsCmd->SetGuidance("kinetic energies of such secondaries are accumulated");
  applyCutsCmd->SetGuidance("to the energy deposition of their mother.");
  applyCutsCmd->SetGuidance(" Note that 'applyCuts' makes sense only for gamma,");
  applyCutsCmd->SetGuidance("e- and e+. If this command is issued for other particle,");
  applyCutsCmd->SetGuidance("a warning message is displayed and the command is");
  applyCutsCmd->SetGuidance("ignored.");
  applyCutsCmd->SetGuidance(" If particle name is 'all', this command affects on");
  applyCutsCmd->SetGuidance("gamma, e- and e+.");
  param = new G4UIparameter("Flag",'s',true);
  param->SetDefaultValue("true");
  applyCutsCmd->SetParameter(param);
  param = new G4UIparameter("Particle",'s',true);
  param->SetDefaultValue("all");
  applyCutsCmd->SetParameter(param);
  applyCutsCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);

  //  /run/particle/dumpCutValues command
  dumpCutValuesCmd = new G4UIcmdWithAString("/run/particle/dumpCutValues",this);
  dumpCutValuesCmd->SetGuidance("Dump a list of production threshold values in range and energy");
  dumpCutValuesCmd->SetGuidance("for all registered material-cuts-couples.");
  dumpCutValuesCmd->SetGuidance("Dumping a list takes place when you issue 'beamOn' and");
  dumpCutValuesCmd->SetGuidance("actual conversion tables from range to energy are available.");
  dumpCutValuesCmd->SetGuidance("If you want a list 'immediately', use '/run/dumpRegion' for threshold");
  dumpCutValuesCmd->SetGuidance("list given in gange only. Also, '/run/dumpCouples' gives you the");
  dumpCutValuesCmd->SetGuidance("current list if you have already issued 'run/beamOn' at least once.");
  dumpCutValuesCmd->SetParameterName("particle",true);
  dumpCutValuesCmd->SetDefaultValue("all");
  dumpCutValuesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

G4UserPhysicsListMessenger::~G4UserPhysicsListMessenger()
{
  delete setPCutCmd; 
  delete setCutCmd; 
  delete setCutRCmd; 
  delete verboseCmd;
  delete dumpListCmd;
  delete addProcManCmd;
  delete buildPTCmd;
  delete storeCmd;  
  delete retrieveCmd;
  delete asciiCmd;
  delete applyCutsCmd;
  delete dumpCutValuesCmd;
  delete theDirectory;
}

void G4UserPhysicsListMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==setCutCmd || command==setPCutCmd ){
    G4double newCut = setCutCmd->GetNewDoubleValue(newValue); 
    thePhysicsList->SetDefaultCutValue(newCut);
    thePhysicsList->SetCutsWithDefault();

  } else if( command==setCutRCmd ){
    //G4Tokenizer next( newValue );
    //G4String rName = G4String(next());
    //G4String cValue = G4String(next())+" "+G4String(next());
    //G4double newCut = setCutCmd->GetNewDoubleValue(cValue); 
    //thePhysicsList->SetCutsForRegion(newCut,rName);
    std::istrstream is(newValue);
    char regName[50];
    G4double cVal;
    char uniName[10];
    is >> regName >> cVal >> uniName;
    G4String regN = regName;
    G4String uniN = uniName;
    thePhysicsList->SetCutsForRegion(cVal*(setCutRCmd->ValueOf(uniN)),regN);

  } else if( command==verboseCmd ) {
    thePhysicsList->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); 

  } else if( command==dumpListCmd ){
    thePhysicsList->DumpList();

  }  else if( command == addProcManCmd ){
    G4ParticleDefinition* particle = (G4ParticleTable::GetParticleTable())->FindParticle(newValue);
    if (particle == 0) return;
    if (particle->GetProcessManager() != 0) return;
    thePhysicsList->AddProcessManager(particle);

  }  else if( command == buildPTCmd ){
    G4ParticleDefinition* particle = (G4ParticleTable::GetParticleTable())->FindParticle(newValue);
    if (particle == 0) return;
    thePhysicsList->PreparePhysicsTable(particle);
    thePhysicsList->BuildPhysicsTable(particle);
    
  } else if ( command == storeCmd ){
    thePhysicsList->StorePhysicsTable(newValue);
  
  } else if( command == retrieveCmd ) {
    if ((newValue == "OFF") || (newValue == "off") ){
      thePhysicsList->ResetPhysicsTableRetrieved();
    } else {
      thePhysicsList->SetPhysicsTableRetrieved(newValue);
    }

  } else if( command == asciiCmd ) {
    if (asciiCmd->GetNewIntValue(newValue) == 0) {
      thePhysicsList->ResetStoredInAscii();
    } else {
      thePhysicsList->SetStoredInAscii();
    }

  } else if( command == applyCutsCmd ) {
    G4Tokenizer next( newValue );

    // check 1st argument
    G4String temp = G4String(next());
    G4bool flag = (temp =="true" || temp=="TRUE");

    // check 2nd argument
    G4String name = G4String(next());

    thePhysicsList->SetApplyCuts(flag, name);
 
  } else if( command == dumpCutValuesCmd ) {
    thePhysicsList->DumpCutValuesTable(1);
  }

} 

G4String G4UserPhysicsListMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  G4String candidates("none");
  G4ParticleTable::G4PTblDicIterator *piter = (G4ParticleTable::GetParticleTable())->GetIterator();
  
  if( command==setCutCmd || command==setPCutCmd ){
    cv = setCutCmd->ConvertToString( thePhysicsList->GetDefaultCutValue(), "mm" );
    
  } else if( command==verboseCmd ){
    cv = verboseCmd->ConvertToString(thePhysicsList->GetVerboseLevel());
    
  }  else if( command== addProcManCmd ){
    // set candidate list
    piter -> reset();
    while( (*piter)() ){
      G4ParticleDefinition *particle = piter->value();
      candidates += " " + particle->GetParticleName();
    }
    addProcManCmd->SetCandidates(candidates);   
    cv = "";
    
  }  else if( command== buildPTCmd ){
    // set candidate list
    piter -> reset();
    while( (*piter)() ){
      G4ParticleDefinition *particle = piter->value();
      candidates += " " + particle->GetParticleName();
    }
    addProcManCmd->SetCandidates(candidates);   
    cv = "";
    
  } else if ( command == storeCmd ){
    cv = thePhysicsList->GetPhysicsTableDirectory();

  }else if( command == retrieveCmd ) {
    if (thePhysicsList->IsPhysicsTableRetrieved()) {
      cv = thePhysicsList->GetPhysicsTableDirectory();
    } else {
      cv = "OFF";
    }

  } else if( command==asciiCmd ){
    if (thePhysicsList->IsStoredInAscii()){
      cv = "1";
    } else {
      cv = "0";
    }

//  } else if( command == applyCutsCmd ) {
//   if (thePhysicsList->GetApplyCuts("gamma")){
//     cv =  "true";
//   } else {
//     cv =  "false";
//   } 
  }
   
  return cv;
}








