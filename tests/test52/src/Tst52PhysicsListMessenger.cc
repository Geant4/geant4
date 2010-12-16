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
// Author: Susanna Guatelli (guatelli@ge.infn.it)
// -------------------------------------------------------------------

#include "Tst52PhysicsListMessenger.hh"
#include "Tst52PhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

Tst52PhysicsListMessenger::Tst52PhysicsListMessenger(Tst52PhysicsList * List)
:Tst52List(List)
{
  EnDir = new G4UIdirectory("/physics/");
  EnDir->SetGuidance("physics commands");

  physicsListCmd = new G4UIcmdWithAString("/physics/addPhysics",this);  
  physicsListCmd->SetGuidance("Add chunks of PhysicsList:photon-standard, photon-epdl, photon-penelope");
  physicsListCmd->SetParameterName("choice",false);
  physicsListCmd->AvailableForStates(G4State_PreInit);  

  
  cutECmd = new G4UIcmdWithADoubleAndUnit("/physics/cutE",this);
  cutECmd->SetGuidance("Set cut values.");
  cutECmd->SetParameterName("range",true);
  cutECmd->SetDefaultValue(1.);
  cutECmd->SetDefaultUnit("mm");
  cutECmd->AvailableForStates(G4State_PreInit);

  cutELowECmd = new G4UIcmdWithADoubleAndUnit("/physics/cutElowe",this);
  cutELowECmd->SetGuidance("Set cut values.");
  cutELowECmd->SetParameterName("energy",true);
  cutELowECmd->SetDefaultValue(1.);
  cutELowECmd->SetDefaultUnit("MeV");
  cutELowECmd->AvailableForStates(G4State_PreInit);

  facRangeCmd = new G4UIcmdWithADouble("/physics/facRange",this);
  facRangeCmd->SetGuidance("Set FacRange.");
  facRangeCmd->SetParameterName("range",true);
  facRangeCmd->AvailableForStates(G4State_PreInit);
}

Tst52PhysicsListMessenger::~Tst52PhysicsListMessenger()
{  
  delete facRangeCmd;
  delete cutELowECmd;
  delete cutECmd;
  delete physicsListCmd;
  delete EnDir;
}
  
void Tst52PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == physicsListCmd)
   { Tst52List->AddPhysicsList(newValue); }

 if (command == cutECmd)
   {
     Tst52List->SetParticleCut(cutECmd->GetNewDoubleValue(newValue)); }
  
 if (command == facRangeCmd)
   {
     G4double facRange = facRangeCmd->GetNewDoubleValue(newValue);
     Tst52List -> SetFacRange(facRange);
   }
}







