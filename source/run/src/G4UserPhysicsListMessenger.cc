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
// $Id: G4UserPhysicsListMessenger.cc 108540 2018-02-16 09:30:43Z gcosmo $
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
//        add dumpOrderingParam command    3 May. 2011 by H.Kurashige
// ------------------------------------------------------------

#include <sstream>

#include "G4UserPhysicsListMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4VUserPhysicsList.hh"
#include "G4PhysicsListHelper.hh"
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
  
  // /run/setCut command
  setCutCmd = new G4UIcmdWithADoubleAndUnit("/run/setCut",this);
  setCutCmd->SetGuidance("Set default cut value ");
  setCutCmd->SetParameterName("cut",false);
  setCutCmd->SetDefaultValue(1.0);
  setCutCmd->SetRange("cut >=0.0");
  setCutCmd->SetDefaultUnit("mm");
  setCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // /run/setCutForAGivenParticle command
  setCutForAGivenParticleCmd = new G4UIcommand("/run/setCutForAGivenParticle",this) ;
  setCutForAGivenParticleCmd->SetGuidance("Set a cut value to a specific particle ") ;
  setCutForAGivenParticleCmd->SetGuidance("Usage: /run/setCutForAGivenParticle  gamma  1. mm") ;
  param = new G4UIparameter("particleName",'s',false) ;
  param->SetParameterCandidates("e- e+ gamma proton");
  setCutForAGivenParticleCmd->SetParameter(param) ;
  param = new G4UIparameter("cut",'d',false) ;
  param->SetDefaultValue("1.") ;
  param->SetParameterRange("cut>=0.0") ;
  setCutForAGivenParticleCmd->SetParameter(param) ;
  param = new G4UIparameter("unit",'s',false) ;
  param->SetDefaultValue("mm") ;
  setCutForAGivenParticleCmd->SetParameter(param) ;
  setCutForAGivenParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // /run/getCutForAGivenParticle command
  getCutForAGivenParticleCmd = new G4UIcmdWithAString("/run/getCutForAGivenParticle",this) ;
  getCutForAGivenParticleCmd->SetGuidance("Get a cut value to a specific particle ") ;
  getCutForAGivenParticleCmd->SetGuidance("Usage: /run/getCutForAGivenParticle  gamma ") ;
  getCutForAGivenParticleCmd->SetParameterName("particleName",false,false) ;
  getCutForAGivenParticleCmd->SetCandidates("e- e+ gamma proton");
  getCutForAGivenParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  // /run/setCutForRegion command
  setCutRCmd = new G4UIcommand("/run/setCutForRegion",this);
  setCutRCmd->SetGuidance("Set cut value for a region");
  param = new G4UIparameter("Region",'s',false);
  setCutRCmd->SetParameter(param);
  param = new G4UIparameter("cut",'d',false);
  param->SetParameterRange("cut >=0.0");
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
  addProcManCmd->SetToBeBroadcasted(false);
  addProcManCmd->SetGuidance("add process manager to specified particle type");
  addProcManCmd->SetParameterName("particleType", true);
  addProcManCmd->SetDefaultValue("");
  addProcManCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

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
  dumpCutValuesCmd->AvailableForStates(G4State_Idle);

  //  /run/particle/dumpCutValues command
  dumpOrdParamCmd = new G4UIcmdWithAnInteger("/run/particle/dumpOrderingParam",this);
  dumpOrdParamCmd->SetGuidance("Dump a list of ordering parameter ");
  dumpOrdParamCmd->SetParameterName("subtype",true);
  dumpOrdParamCmd->SetDefaultValue(-1);
  dumpOrdParamCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);
}

G4UserPhysicsListMessenger::~G4UserPhysicsListMessenger()
{
  delete setCutCmd; 
  delete setCutRCmd;
  delete setCutForAGivenParticleCmd;
  delete getCutForAGivenParticleCmd;
  delete verboseCmd;
  delete dumpListCmd;
  delete addProcManCmd;
  delete buildPTCmd;
  delete storeCmd;  
  delete retrieveCmd;
  delete asciiCmd;
  delete applyCutsCmd;
  delete dumpCutValuesCmd;
  delete dumpOrdParamCmd;
  delete theDirectory;
}

void G4UserPhysicsListMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  G4ExceptionDescription ed;
  if( command==setCutCmd ){
    G4double newCut = setCutCmd->GetNewDoubleValue(newValue); 
    thePhysicsList->SetDefaultCutValue(newCut);
    thePhysicsList->SetCuts();

  } else if( command==setCutForAGivenParticleCmd ){
    G4String particleName, unit ; G4double cut ;
    std::istringstream str (newValue) ;
    str >> particleName >> cut >> unit ;
    thePhysicsList->SetCutValue(cut*G4UIcommand::ValueOf(unit), particleName) ; 

  } else if( command==getCutForAGivenParticleCmd ){
    G4cout << thePhysicsList->GetCutValue(newValue)/mm <<"[mm]" << G4endl ;

  } else if( command==setCutRCmd ){
    std::istringstream is(newValue);
    G4String regName;
    G4String uniName;
    G4double cVal = -1.0;
    is >> regName >> cVal >> uniName;
    if (is.fail()) {
      ed << "illegal arguments : " << newValue;
      command->CommandFailed(ed);
      return;
    }
    thePhysicsList->SetCutsForRegion(cVal*(setCutRCmd->ValueOf(uniName)),regName);

  } else if( command==verboseCmd ) {
    thePhysicsList->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); 

  } else if( command==dumpListCmd ){
    thePhysicsList->DumpList();

  } else if( command==dumpOrdParamCmd ){
    G4int stype = dumpOrdParamCmd->GetNewIntValue(newValue);
    G4PhysicsListHelper::GetPhysicsListHelper()->DumpOrdingParameterTable(stype);

  }  else if( command == addProcManCmd ){
    G4ParticleDefinition* particle = (G4ParticleTable::GetParticleTable())->FindParticle(newValue);
    if (particle == 0)
    {
      ed << " Particle is not found : " << newValue;
      command->CommandFailed(ed);
      return;
    }
    else if (particle->GetProcessManager() != 0)
    {
      ed << " Particle is not initialized : " << newValue;
      command->CommandFailed(ed);
      return;
    }
    thePhysicsList->AddProcessManager(particle);

  }  else if( command == buildPTCmd ){
    G4ParticleDefinition* particle = (G4ParticleTable::GetParticleTable())->FindParticle(newValue);
    if (particle == 0)
    {
      ed << " Particle is not found : " << newValue;
      command->CommandFailed(ed);
      return;
    }
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
  
  if( command==setCutCmd ) {
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
