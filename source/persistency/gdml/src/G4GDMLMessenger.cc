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
// $Id: G4GDMLMessenger.cc 105114 2017-07-13 09:49:54Z gcosmo $
//
//
// class G4GDMLMessenger Implementation
//
// -------------------------------------------------------------------------

#include "G4GDMLMessenger.hh"
#include "G4GDMLParser.hh"

#include "globals.hh"
#include "G4RunManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"

G4GDMLMessenger::G4GDMLMessenger(G4GDMLParser* myPars)
  : myParser(myPars), topvol(0), pFlag(true)
{ 
  persistencyDir = new G4UIdirectory("/persistency/");
  persistencyDir->SetGuidance("UI commands specific to persistency.");
  
  gdmlDir = new G4UIdirectory("/persistency/gdml/");
  gdmlDir->SetGuidance("GDML parser and writer.");
  
  ReaderCmd = new G4UIcmdWithAString("/persistency/gdml/read",this);
  ReaderCmd->SetGuidance("Read GDML file.");
  ReaderCmd->SetParameterName("filename",false);
  ReaderCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TopVolCmd = new G4UIcmdWithAString("/persistency/gdml/topvol",this);
  TopVolCmd->SetGuidance("Set the top volume for writing the GDML file.");
  TopVolCmd->SetParameterName("topvol",false);

  WriterCmd = new G4UIcmdWithAString("/persistency/gdml/write",this);
  WriterCmd->SetGuidance("Write GDML file.");
  WriterCmd->SetParameterName("filename",false);
  WriterCmd->AvailableForStates(G4State_Idle);

  StripCmd = new G4UIcmdWithABool("/persistency/gdml/strip_pointers",this);
  StripCmd->SetGuidance("Enable/disable stripping of pointers on names");
  StripCmd->SetGuidance("when reading a GDML file.");
  StripCmd->SetParameterName("strip_pointers",true);
  StripCmd->SetDefaultValue(true);
  StripCmd->AvailableForStates(G4State_Idle);

  AppendCmd = new G4UIcmdWithABool("/persistency/gdml/add_pointers",this);
  AppendCmd->SetGuidance("Enable/disable appending of pointers to names");
  AppendCmd->SetGuidance("when writing a GDML file.");
  AppendCmd->SetParameterName("append_pointers",true);
  AppendCmd->SetDefaultValue(true);
  AppendCmd->AvailableForStates(G4State_Idle);

  RegionCmd = new G4UIcmdWithABool("/persistency/gdml/export_regions",this);
  RegionCmd->SetGuidance("Enable export of geometrical regions");
  RegionCmd->SetGuidance("for storing production cuts.");
  RegionCmd->SetParameterName("export_regions",false);
  RegionCmd->SetDefaultValue(false);
  RegionCmd->AvailableForStates(G4State_Idle);

  EcutsCmd = new G4UIcmdWithABool("/persistency/gdml/export_Ecuts",this);
  EcutsCmd->SetGuidance("Enable export of energy cuts associated");
  EcutsCmd->SetGuidance("to logical volumes.");
  EcutsCmd->SetGuidance("NOTE: may increase considerably the size of the");
  EcutsCmd->SetGuidance("      GDML file!  Information is anyhow not used");
  EcutsCmd->SetGuidance("      for import.");
  EcutsCmd->SetParameterName("export_Ecuts",false);
  EcutsCmd->SetDefaultValue(false);
  EcutsCmd->AvailableForStates(G4State_Idle);

  SDCmd = new G4UIcmdWithABool("/persistency/gdml/export_SD",this);
  SDCmd->SetGuidance("Enable export of SD associated");
  SDCmd->SetGuidance("to logical volumes.");
  SDCmd->SetParameterName("export_SD",false);
  SDCmd->SetDefaultValue(false);
  SDCmd->AvailableForStates(G4State_Idle);

  ClearCmd = new G4UIcmdWithoutParameter("/persistency/gdml/clear",this);
  ClearCmd->SetGuidance("Clear geometry (before reading a new one from GDML).");
  ClearCmd->AvailableForStates(G4State_Idle);
}

G4GDMLMessenger::~G4GDMLMessenger()
{
  delete ReaderCmd;
  delete WriterCmd;
  delete ClearCmd;
  delete TopVolCmd;
  delete RegionCmd;
  delete EcutsCmd;
  delete SDCmd;
  delete persistencyDir;
  delete gdmlDir;
  delete StripCmd;
  delete AppendCmd;
}

void G4GDMLMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == StripCmd )
  {
    G4bool mode = StripCmd->GetNewBoolValue(newValue);
    myParser->SetStripFlag(mode);
  }

  if( command == AppendCmd )
  {
    pFlag = AppendCmd->GetNewBoolValue(newValue);
    myParser->SetAddPointerToName(pFlag);
  }

  if( command == ReaderCmd )
  { 
    G4GeometryManager::GetInstance()->OpenGeometry();
    myParser->Read(newValue);
    G4RunManager::GetRunManager()->DefineWorldVolume(myParser->GetWorldVolume());
  }

  if( command == RegionCmd )
  {
    G4bool mode = RegionCmd->GetNewBoolValue(newValue);
    myParser->SetRegionExport(mode);
  }

  if( command == EcutsCmd )
  {
    G4bool mode = EcutsCmd->GetNewBoolValue(newValue);
    myParser->SetEnergyCutsExport(mode);
  }
   
  if( command == SDCmd )
  {
    G4bool mode = SDCmd->GetNewBoolValue(newValue);
    myParser->SetSDExport(mode);
  }
   
  if( command == TopVolCmd )
  {
    topvol = G4LogicalVolumeStore::GetInstance()->GetVolume(newValue);
  }
   
  if( command == WriterCmd )
  {
    myParser->Write(newValue, topvol, pFlag);
  }  

  if( command == ClearCmd )
  { 
    myParser->Clear();
    G4RunManager::GetRunManager()->ReinitializeGeometry(true);
  }  
}
