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
#include "G4HadronicEPTestMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4HadronicException.hh"


G4HadronicEPTestMessenger::G4HadronicEPTestMessenger(G4HadronicProcessStore* theStore)
 :theProcessStore(theStore)
{
  // Main directory for control of the e/p test
  old_heptstDirectory = new G4UIdirectory("/heptst/");                                // To be removed in G4 11.0
  old_heptstDirectory->SetGuidance("Controls for the hadronic energy/momentum test");

  heptstDirectory = new G4UIdirectory("/process/had/heptst/");
  heptstDirectory->SetGuidance("Controls for the hadronic energy/momentum test");
  
  // Command to set level of detail reported upon e/p non-conservation
  old_reportLvlCmd = new G4UIcmdWithAnInteger("/heptst/reportLevel",this);            // To be removed in G4 11.0
  old_reportLvlCmd->SetGuidance("Set level of detail reported upon E/p non-conservation");
  old_reportLvlCmd->SetGuidance(" 0 - (default) no reporting "); 
  old_reportLvlCmd->SetGuidance(" 1 - report only when E/p not conserved "); 
  old_reportLvlCmd->SetGuidance(" 2 - report regardless of E/p conservation "); 
  old_reportLvlCmd->SetGuidance(" 3 - report only when E/p not conserved, with names, limits "); 
  old_reportLvlCmd->SetGuidance(" 4 - report regardless of E/p conservation, with names, limits "); 
  old_reportLvlCmd->SetParameterName("OldReportLevel",true);
  old_reportLvlCmd->SetDefaultValue(0);
  old_reportLvlCmd->SetRange("OldReportLevel >= 0 && OldReportLevel < 5");
  
  reportLvlCmd = new G4UIcmdWithAnInteger("/process/had/heptst/reportLevel",this);
  reportLvlCmd->SetGuidance("Set level of detail reported upon E/p non-conservation");
  reportLvlCmd->SetGuidance(" 0 - (default) no reporting "); 
  reportLvlCmd->SetGuidance(" 1 - report only when E/p not conserved "); 
  reportLvlCmd->SetGuidance(" 2 - report regardless of E/p conservation "); 
  reportLvlCmd->SetGuidance(" 3 - report only when E/p not conserved, with names, limits "); 
  reportLvlCmd->SetGuidance(" 4 - report regardless of E/p conservation, with names, limits "); 
  reportLvlCmd->SetParameterName("ReportLevel",true);
  reportLvlCmd->SetDefaultValue(0);
  reportLvlCmd->SetRange("ReportLevel >= 0 && ReportLevel < 5");
  
  // Set the relative energy non-conservation level for the process
  old_procRelLvlCmd = new G4UIcmdWithADouble("/heptst/processRelLevel",this);         // To be removed in G4 11.0
  old_procRelLvlCmd->SetGuidance("Set relative level of allowed energy non-conservation");
  old_procRelLvlCmd->SetParameterName("OlProcessRelLevel",true);
  old_procRelLvlCmd->SetDefaultValue(-1.0);

  procRelLvlCmd = new G4UIcmdWithADouble("/process/had/heptst/processRelLevel",this);
  procRelLvlCmd->SetGuidance("Set relative level of allowed energy non-conservation");
  procRelLvlCmd->SetParameterName("ProcessRelLevel",true);
  procRelLvlCmd->SetDefaultValue(-1.0);

  // Set the absolute energy non-conservation level for the process
  old_procAbsLvlCmd = new G4UIcmdWithADoubleAndUnit("/heptst/processAbsLevel",this);  // To be removed in G4 11.0
  old_procAbsLvlCmd->SetGuidance("Set absolute energy level (with unit) of allowed energy non-conservation");
  old_procAbsLvlCmd->SetParameterName("OldProcessAbsLevel",true);
  old_procAbsLvlCmd->SetDefaultValue(-1.0);
  old_procAbsLvlCmd->SetUnitCategory("Energy");

  procAbsLvlCmd = new G4UIcmdWithADoubleAndUnit("/process/had/heptst/processAbsLevel",this);
  procAbsLvlCmd->SetGuidance("Set absolute energy level (with unit) of allowed energy non-conservation");
  procAbsLvlCmd->SetParameterName("ProcessAbsLevel",true);
  procAbsLvlCmd->SetDefaultValue(-1.0);
  procAbsLvlCmd->SetUnitCategory("Energy");
}


G4HadronicEPTestMessenger::~G4HadronicEPTestMessenger()
{
  delete old_heptstDirectory;  // To be removed in G4 11.0
  delete     heptstDirectory;
  delete old_reportLvlCmd;     // To be removed in G4 11.0
  delete     reportLvlCmd;
  delete old_procRelLvlCmd;    // To be removed in G4 11.0
  delete     procRelLvlCmd;
  delete old_procAbsLvlCmd;    // To be removed in G4 11.0
  delete     procAbsLvlCmd; 
}


void G4HadronicEPTestMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  // Old commands to be removed in G4 11.0
  if ( command == old_reportLvlCmd ) {
    theProcessStore->SetEpReportLevel( old_reportLvlCmd->GetNewIntValue( newValue ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/heptst/reportLevel in the next major release, Geant4 version 11.0";
    G4Exception( "G4HadronicEPTestMessenger", "hadEPTestMessenger001", JustWarning, ed );    
  } else if ( command == old_procRelLvlCmd ) {
    theProcessStore->SetProcessRelLevel( old_procRelLvlCmd->GetNewDoubleValue( newValue ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/heptst/processRelLevel in the next major release, Geant4 version 11.0";
    G4Exception( "G4HadronicEPTestMessenger", "hadEPTestMessenger002", JustWarning, ed );
  } else if ( command == old_procAbsLvlCmd ) {
    theProcessStore->SetProcessAbsLevel( old_procAbsLvlCmd->GetNewDoubleValue( newValue ) );
    G4ExceptionDescription ed;
    ed << "This command is valid but deprecated and will be replaced with the command:\n"
       << "/process/had/heptst/processAbsLevel in the next major release, Geant4 version 11.0";
    G4Exception( "G4HadronicEPTestMessenger", "hadEPTestMessenger003", JustWarning, ed );
  }
  
  // New commands
  if ( command == reportLvlCmd ) {
    theProcessStore->SetEpReportLevel( reportLvlCmd->GetNewIntValue( newValue ) );
  } else if ( command == procRelLvlCmd ) {
    theProcessStore->SetProcessRelLevel( procRelLvlCmd->GetNewDoubleValue( newValue ) );
  } else if ( command == procAbsLvlCmd ) {
    theProcessStore->SetProcessAbsLevel( procAbsLvlCmd->GetNewDoubleValue( newValue ) );
  }
}
