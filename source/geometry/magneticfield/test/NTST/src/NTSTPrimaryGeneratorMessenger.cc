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
// $Id: NTSTPrimaryGeneratorMessenger.cc,v 1.2 2003-12-09 15:35:41 gunter Exp $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "NTSTPrimaryGeneratorMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTPrimaryGeneratorMessenger::NTSTPrimaryGeneratorMessenger()
  : ChooseCmd(0), PrintCmd(0), generator(0), print(0)
{ // set the available commands
  GenDir = new G4UIdirectory("/gen/");
  GenDir->SetGuidance("NTST generator control");
  ChooseCmd = new G4UIcmdWithAString("/gen/choose",this);
  ChooseCmd->SetGuidance("Choose the generator to use.");
  ChooseCmd->SetDefaultValue("gun");
  G4String Names("gun evt");
  //  ChooseCmd->SetCandidates(Names);
  //  ChooseCmd->AvailableForStates(PreInit,Idle);
  SetNewValue(ChooseCmd, "gun");

  // (en,dis)able printing
  PrintCmd = new G4UIcmdWithABool("/gen/print", this);
  G4cout << "### Instantiated /gen/print cmd" << G4endl;
  PrintCmd->SetGuidance("print events (true) or not (false).");
  G4cout << "### SetGuidance for /gen/print" << G4endl;
  PrintCmd->SetDefaultValue(false);
  G4cout << "### Set default value for /gen/print cmd" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTPrimaryGeneratorMessenger::~NTSTPrimaryGeneratorMessenger()
{
  delete ChooseCmd;
  delete PrintCmd;
  delete generator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "NTSTGunGenerator.hh"
#include "NTSTBabarEvtReadGenerator.hh"

void 
NTSTPrimaryGeneratorMessenger::
SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == ChooseCmd ) {
    if (newValue == "gun") {
      Choice = GUN;
      delete generator;
      generator = new NTSTGunGenerator();
      G4cout << "### instantiated " << newValue << " generator" << G4endl;
      Name = newValue;
    } else if (newValue == "evt") {
      Choice = EVT;
      delete generator;
      generator = new NTSTBabarEvtReadGenerator("babarevt.out");
      G4cout << "### instantiated " << newValue << " generator" << G4endl;
      Name = newValue;
    } else {
      G4cerr << "!!! Illegal choice of generator: \"" << newValue 
	     << "\"" << G4endl;
      G4cerr << "!!! Valid choices are " << GetNames() << G4endl;
      G4cerr << "!!! Old generator: \"" << GetName() << "\" retained" 
	     << G4endl;
    }
  } else if ( command == PrintCmd ) {
    if (newValue == "0" || newValue == "false") DisablePrinting();
    else if (newValue == "1" || newValue == "true") EnablePrinting();
    else G4cerr << "PrintCmd \"" << newValue << "\" not understood" << G4endl;
  }
}




