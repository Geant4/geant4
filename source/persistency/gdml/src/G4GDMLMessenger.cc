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
// $Id: G4GDMLMessenger.cc,v 1.3 2010-11-02 11:44:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4UIcmdWithoutParameter.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"

G4GDMLMessenger::G4GDMLMessenger(G4GDMLParser* myPars)
  : myParser(myPars), topvol(0)
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
  delete persistencyDir;
  delete gdmlDir;
}

void G4GDMLMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == ReaderCmd )
  { 
    G4GeometryManager::GetInstance()->OpenGeometry();
    myParser->Read(newValue);
    G4RunManager::GetRunManager()->DefineWorldVolume(myParser->GetWorldVolume());
  }

  if( command == TopVolCmd )
  {
    topvol = G4LogicalVolumeStore::GetInstance()->GetVolume(newValue);
  }
   
  if( command == WriterCmd )
  {
    myParser->Write(newValue, topvol);
  }  

  if( command == ClearCmd )
  { 
    myParser->Clear();
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();
  }  
}
