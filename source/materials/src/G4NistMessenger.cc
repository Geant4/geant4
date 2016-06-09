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
// $Id: G4NistMessenger.cc,v 1.4 2007/05/02 10:48:52 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//
// File name:     G4NistMessenger
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.12.2004
//
// Modifications:
//
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4NistMessenger.hh"

#include "G4NistManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4NistMessenger::G4NistMessenger(G4NistManager* man)
:manager(man)
{
  matDir = new G4UIdirectory("/material/");
  matDir->SetGuidance("Commands for materials");

  verCmd = new G4UIcmdWithAnInteger("/material/verbose",this);
  verCmd->SetGuidance("Set verbose level.");
  
  nistDir = new G4UIdirectory("/material/nist/");
  nistDir->SetGuidance("Commands for the nist dataBase");
    
  prtElmCmd = new G4UIcmdWithAString("/material/nist/printElement",this);
  prtElmCmd->SetGuidance("print element(s) in dataBase.");
  prtElmCmd->SetGuidance("symbol = element.");
  prtElmCmd->SetGuidance("all    = all elements.");
  prtElmCmd->SetGuidance("verbose>1 : list associated isotopes.");  
  prtElmCmd->SetParameterName("symbol", true);
  prtElmCmd->SetDefaultValue("all");
  
  przElmCmd = new G4UIcmdWithAnInteger("/material/nist/printElementZ",this);
  przElmCmd->SetGuidance("print element Z in dataBase.");
  przElmCmd->SetGuidance("verbose>1 : list associated isotopes.");
   
  lisMatCmd = new G4UIcmdWithAString("/material/nist/listMaterials",this);
  lisMatCmd->SetGuidance("list materials in dataBase.");
  lisMatCmd->SetGuidance("simple = simple nist materials.");
  lisMatCmd->SetGuidance("compound = compound nist materials.");
  lisMatCmd->SetGuidance("hep = hep materials.");
  lisMatCmd->SetGuidance("all = all materials.");
  lisMatCmd->SetParameterName("list", true);
  lisMatCmd->SetCandidates("simple compound hep all");
  lisMatCmd->SetDefaultValue("all");
  
  g4Dir = new G4UIdirectory("/material/g4/");
  g4Dir->SetGuidance("Commands for G4MaterialTable");
  
  g4ElmCmd = new G4UIcmdWithAString("/material/g4/printElement",this);
  g4ElmCmd->SetGuidance("print Element in G4ElementTable.");
    
  g4MatCmd = new G4UIcmdWithAString("/material/g4/printMaterial",this);
  g4MatCmd->SetGuidance("print Material in G4MaterialTable.");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4NistMessenger::~G4NistMessenger()
{
  delete verCmd;  
  delete prtElmCmd;
  delete przElmCmd;
  delete lisMatCmd;
  delete nistDir;
  
  delete g4ElmCmd;   
  delete g4MatCmd;
  delete g4Dir;
  delete matDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == verCmd) 
   { manager->SetVerbose(verCmd->GetNewIntValue(newValue));}
    
  if (command == prtElmCmd)
   { manager->PrintElement(newValue);}
    
  if (command == przElmCmd)
   { manager->PrintElement(przElmCmd->GetNewIntValue(newValue));}
   
  if (command == lisMatCmd) 
   { manager->ListMaterials(newValue);}

  if (command == g4ElmCmd)
   { manager->PrintG4Element(newValue);}
   
  if (command == g4MatCmd)
   { manager->PrintG4Material(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
