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
// $Id: MedLinacPhysicsListMessenger.cc,v 1.2 2006/06/29 16:04:41 gunter Exp $
//
//  Code developed by: M. Piergentili

#include "MedLinacPhysicsListMessenger.hh"

#include "MedLinacPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//*********************************************************************

MedLinacPhysicsListMessenger::MedLinacPhysicsListMessenger(
                                           MedLinacPhysicsList* MedLinacCut)
:pMedLinacPhysicsList(MedLinacCut)
{ 

  G4cout <<"==================PhysicsListMessenger  "<<G4endl;
  PhysicsDir = new G4UIdirectory("/PhysicsList/");
  PhysicsDir->SetGuidance("physics parameters");

  CutCmd = new G4UIcmdWithADoubleAndUnit("/PhysicsList/cut",this);
  CutCmd->SetGuidance("Set cut length (mm)");
  CutCmd->SetParameterName("defaultCut",false);
  CutCmd->SetDefaultUnit( "mm" );
  CutCmd->SetUnitCategory("Length");
  CutCmd->AvailableForStates(G4State_PreInit);

  //UpdateCmd = new G4UIcmdWithoutParameter("/Jaws/update",this);
  //UpdateCmd->SetGuidance("Update geometry.");
  //UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  //UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  //UpdateCmd->AvailableForStates(G4State_Idle);


}

//****************************************************************************
MedLinacPhysicsListMessenger::~MedLinacPhysicsListMessenger()
{

  delete CutCmd;
  //delete UpdateCmd;

  delete PhysicsDir;
}

//****************************************************************************

void MedLinacPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == CutCmd )
   { pMedLinacPhysicsList->SetCut(CutCmd->GetNewDoubleValue(newValue));}

  //if( command == UpdateCmd )
  // { MedLinacDetector->UpdateGeometry(); }

}

//****************************************************************************
