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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraRunActionMessenger.cc
//    ****************************************************
//
//    Messenger Class for UltraRunAction
//    Allows to set the run ID
//
#include "UltraRunActionMessenger.hh"
#include "UltraRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include <sstream>

UltraRunActionMessenger::UltraRunActionMessenger(UltraRunAction* aRunAction)
:theRunAction(aRunAction)
{
  runDirectory = new G4UIdirectory("/mysetrun/");
  runDirectory->SetGuidance("My set commands.");

  runIDCmd = new G4UIcommand("/mysetrun/SetRunID",this);
  runIDCmd->SetGuidance("Set run ID");
  runIDCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* p1 = new G4UIparameter("runID",'i',true);
  p1->SetDefaultValue(1000);
  runIDCmd->SetParameter(p1);
}

UltraRunActionMessenger::~UltraRunActionMessenger()
{
  delete runIDCmd;
  delete runDirectory;
}

void UltraRunActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  const char* nv = (const char*)newValue;
  if( command==runIDCmd )
  {
    G4int id;
    std::istringstream is(nv);
    is >> id;

    theRunAction->MySetRunID(id);  
  }
}



