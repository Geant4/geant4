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
// $Id: MedLinacMLCMessenger.cc,v 1.3 2006/06/29 16:04:29 gunter Exp $
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
