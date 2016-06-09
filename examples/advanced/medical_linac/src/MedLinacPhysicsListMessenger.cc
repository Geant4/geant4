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
// $Id: MedLinacPhysicsListMessenger.cc,v 1.1 2005/07/03 23:27:37 mpiergen Exp $
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
