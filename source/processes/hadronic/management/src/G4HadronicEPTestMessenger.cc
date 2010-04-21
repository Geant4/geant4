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


G4HadronicEPTestMessenger::G4HadronicEPTestMessenger(G4HadronicProcessStore* theStore)
 :theProcessStore(theStore)
{
  // Main directory for control of the e/p test
  heptstDirectory = new G4UIdirectory("/heptst/");
  heptstDirectory->SetGuidance("Controls for the hadronic energy/momentum test");

  // Command to set level of detail reported upon e/p non-conservation
  reportLvlCmd = new G4UIcmdWithAnInteger("/heptst/reportLevel",this);
  reportLvlCmd->SetGuidance("Set level of detail reported upon E/p non-conservation");
  reportLvlCmd->SetGuidance(" 0 - (default) no reporting "); 
  reportLvlCmd->SetGuidance(" 1 - report only when E/p not conserved "); 
  reportLvlCmd->SetGuidance(" 2 - report regardless of E/p conservation "); 
  reportLvlCmd->SetGuidance(" 3 - report only when E/p not conserved, with names, limits "); 
  reportLvlCmd->SetGuidance(" 4 - report regardless of E/p conservation, with names, limits "); 
  reportLvlCmd->SetParameterName("ReportLevel",true);
  reportLvlCmd->SetDefaultValue(0);
  reportLvlCmd->SetRange("ReportLevel >= 0 && ReportLevel < 5");

  // Set the absolute energy non-conservation level for the process
  procAbsLvlCmd = new G4UIcmdWithADoubleAndUnit("/heptst/processAbsLevel",this);
  procAbsLvlCmd->SetGuidance("Set absolute energy level (with unit) of allowed energy non-conservation");
  procAbsLvlCmd->SetParameterName("ProcessAbsLevel",true);
  procAbsLvlCmd->SetDefaultValue(-1.0);
  procAbsLvlCmd->SetUnitCategory("Energy");

  // Set the relative energy non-conservation level for the process
  procRelLvlCmd = new G4UIcmdWithADouble("/heptst/processRelLevel",this);
  procRelLvlCmd->SetGuidance("Set relative level of allowed energy non-conservation");
  procRelLvlCmd->SetParameterName("ProcessRelLevel",true);
  procRelLvlCmd->SetDefaultValue(-1.0);
}


G4HadronicEPTestMessenger::~G4HadronicEPTestMessenger()
{
  delete heptstDirectory;
  delete reportLvlCmd;
  delete procAbsLvlCmd;
  delete procRelLvlCmd;
}


void G4HadronicEPTestMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command==reportLvlCmd) {
    theProcessStore->SetEpReportLevel(reportLvlCmd->GetNewIntValue(newValue) );

  } else if(command==procRelLvlCmd) {
    //    G4double absval = theHadronicProcess->GetEnergyMomentumCheckLevels().second;
    theProcessStore->SetProcessRelLevel(procRelLvlCmd->GetNewDoubleValue(newValue) );
  } else if(command==procAbsLvlCmd) {
    //    G4double relval = theHadronicProcess->GetEnergyMomentumCheckLevels().first;
    theProcessStore->SetProcessAbsLevel(procAbsLvlCmd->GetNewDoubleValue(newValue) );
  }
}
