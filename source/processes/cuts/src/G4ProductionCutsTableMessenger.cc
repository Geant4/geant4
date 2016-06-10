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
// $Id: G4ProductionCutsTableMessenger.cc 70369 2013-05-29 14:59:24Z gcosmo $
//
// 
//---------------------------------------------------------------
//
//  G4ProductionCutsTableMessenger.cc
// ------------------------------------------------------------
//      History
//        first version                   02 Mar. 2008 by H.Kurashige 
//

#include "G4ProductionCutsTableMessenger.hh"
#include "G4ProductionCutsTable.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#include "G4Tokenizer.hh"           

#include <sstream>

G4ProductionCutsTableMessenger::G4ProductionCutsTableMessenger(   G4ProductionCutsTable* pTable)
  :theCutsTable(pTable)
{
  // /cuts/   directory
  theDirectory = new G4UIdirectory("/cuts/");
  theDirectory->SetGuidance("Commands for G4VUserPhysicsList.");

  // /cuts/verbose command
  verboseCmd = new G4UIcmdWithAnInteger("/cuts/verbose",this);
  verboseCmd->SetGuidance("Set the Verbose level of G4ProductionCutsTable.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display warning messages");
  verboseCmd->SetGuidance(" 2 : Display more info");
  verboseCmd->SetGuidance(" 2 : Display debug info");
  verboseCmd->SetParameterName("level",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("level >=0 && level <=3");

  // /cuts/setLowEdge command
  setLowEdgeCmd = new G4UIcmdWithADoubleAndUnit("/cuts/setLowEdge",this);
  setLowEdgeCmd->SetGuidance("Set low edge energy value ");
  setLowEdgeCmd->SetParameterName("edge",false);
  setLowEdgeCmd->SetDefaultValue(0.99);
  setLowEdgeCmd->SetRange("edge >0.0");
  setLowEdgeCmd->SetDefaultUnit("keV");
  setLowEdgeCmd->AvailableForStates(G4State_PreInit);

  // /cuts/setHighEdge command
  setHighEdgeCmd = new G4UIcmdWithADoubleAndUnit("/cuts/setHighEdge",this);
  setHighEdgeCmd->SetGuidance("Set high edge energy value ");
  setHighEdgeCmd->SetParameterName("edge",false);
  setHighEdgeCmd->SetDefaultValue(100.0);
  setHighEdgeCmd->SetRange("edge >0.0");
  setHighEdgeCmd->SetDefaultUnit("TeV");
  setHighEdgeCmd->AvailableForStates(G4State_PreInit);
 
  // /cuts/setMaxCutEnergy command
  setMaxEnergyCutCmd = new G4UIcmdWithADoubleAndUnit("/cuts/setMaxCutEnergy",this);
  setMaxEnergyCutCmd->SetGuidance("Set maximum of cut energy value ");
  setMaxEnergyCutCmd->SetParameterName("cut",false);
  setMaxEnergyCutCmd->SetDefaultValue(10.0);
  setMaxEnergyCutCmd->SetRange("cut >0.0");
  setMaxEnergyCutCmd->SetDefaultUnit("GeV");
  setMaxEnergyCutCmd->AvailableForStates(G4State_PreInit);
 
 // /cuts/dump command
  dumpCmd = new G4UIcmdWithoutParameter("/cuts/dump",this);
  dumpCmd->SetGuidance("Dump cuplues in ProductuinCutsTable. ");

}

G4ProductionCutsTableMessenger::~G4ProductionCutsTableMessenger()
{
  delete dumpCmd;
  delete setMaxEnergyCutCmd;
  delete setHighEdgeCmd;
  delete setLowEdgeCmd;
  delete verboseCmd;
  delete theDirectory;
}

void G4ProductionCutsTableMessenger::SetNewValue(G4UIcommand * command,
						 G4String newValue)
{
  if( command==verboseCmd ) {
    theCutsTable->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); 

  } else if( command==dumpCmd ){
    theCutsTable-> DumpCouples();

  } else if( command==setLowEdgeCmd ){
    G4double lowEdge = setLowEdgeCmd->GetNewDoubleValue(newValue); 
    G4double highEdge = theCutsTable->GetHighEdgeEnergy();
    theCutsTable->SetEnergyRange(lowEdge, highEdge);
    
  } else if( command==setHighEdgeCmd ){
    G4double highEdge = setHighEdgeCmd->GetNewDoubleValue(newValue); 
    G4double lowEdge = theCutsTable->GetLowEdgeEnergy();
    theCutsTable->SetEnergyRange(lowEdge, highEdge);
    
 } else if( command==setMaxEnergyCutCmd ){
    G4double cut = setHighEdgeCmd->GetNewDoubleValue(newValue); 
    theCutsTable->SetMaxEnergyCut(cut);
    
  }
}

G4String G4ProductionCutsTableMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==verboseCmd ){
   cv = verboseCmd->ConvertToString(theCutsTable->GetVerboseLevel());

 } else if( command==setLowEdgeCmd ){
    G4double lowEdge = theCutsTable->GetLowEdgeEnergy();
    cv = setLowEdgeCmd->ConvertToString( lowEdge, "keV" );

 } else if( command==setHighEdgeCmd ){
    G4double highEdge = theCutsTable->GetHighEdgeEnergy();
    cv = setHighEdgeCmd->ConvertToString( highEdge, "TeV" );

 } else if( command==setMaxEnergyCutCmd ){
    G4double cut = theCutsTable->GetMaxEnergyCut();
    cv = setMaxEnergyCutCmd->ConvertToString( cut, "GeV" );
 }


  return cv;

}

