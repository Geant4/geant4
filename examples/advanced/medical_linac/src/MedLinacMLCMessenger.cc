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
// $Id: MedLinacMLCMessenger.cc,v 1.2 2004/11/25 16:28:39 mpiergen Exp $
//
//  Code developed by: M. Piergentili

#include "MedLinacMLCMessenger.hh"

#include "MedLinacMLCDecorator.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//*********************************************************************

MedLinacMLCMessenger::MedLinacMLCMessenger(
                                           MedLinacMLCDecorator* MedLinacMLC)
:MedLinacMLCDeco(MedLinacMLC)
{ 

  MedLinacMLCDir = new G4UIdirectory("/MLC/");
  MedLinacMLCDir->SetGuidance("A1 position");
  
  
  leafNameCmd = new G4UIcmdWithAString("/MLC/leaf_selection",this);
  leafNameCmd->SetGuidance("Select the leaf ");
  leafNameCmd->SetParameterName("leafName",false);
  leafNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  positionCmd = new G4UIcmdWithADoubleAndUnit("/MLC/position",this);
  positionCmd->SetGuidance("Set position ");
  positionCmd->SetParameterName("Pos_y",false);
  positionCmd->SetRange("Pos_y>=0. && Pos_y<15.0");
  positionCmd->SetDefaultUnit( "cm" );
  positionCmd->SetUnitCategory("Length");
  positionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}
//****************************************************************************
MedLinacMLCMessenger::~MedLinacMLCMessenger()
{
  delete leafNameCmd;
  delete positionCmd;
  delete MedLinacMLCDir;

}
//****************************************************************************
void MedLinacMLCMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == leafNameCmd )
   MedLinacMLCDeco->SetLeafName(newValue);
  if( command == positionCmd )
   MedLinacMLCDeco->SetPos_y(positionCmd->GetNewDoubleValue(newValue));
 }


//****************************************************************************
