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
#include <strstream>

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
    std::istrstream is((char*)nv);
    is >> id;

    theRunAction->MySetRunID(id);  
  }
}



