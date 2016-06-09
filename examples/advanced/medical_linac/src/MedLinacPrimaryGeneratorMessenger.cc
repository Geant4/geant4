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
// $Id: MedLinacPrimaryGeneratorMessenger.cc,v 1.3 2004/05/14 18:25:40 mpiergen Exp $
//
//
// Code developed by: M. Piergentili
//
#include "MedLinacPrimaryGeneratorMessenger.hh"

#include "MedLinacPrimaryGeneratorAction.hh"
#include "MedLinacRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

MedLinacPrimaryGeneratorMessenger::MedLinacPrimaryGeneratorMessenger(MedLinacPrimaryGeneratorAction* MedLinacGun)
  :MedLinacAction(MedLinacGun)

{  
  EnergyCmd = new G4UIcmdWithADoubleAndUnit("/energy",this);
  EnergyCmd->SetGuidance("Select the energy (and unit) of the primary particles");
  EnergyCmd->SetGuidance("Type particle Energy");
  EnergyCmd->SetParameterName("energy",true,true);
  EnergyCmd->SetDefaultValue(1.);
  EnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SourceTypeCmd = new G4UIcmdWithADoubleAndUnit("/sourceType",this);
  SourceTypeCmd->SetGuidance("Select the standard deviation of the energy.");
  SourceTypeCmd->SetParameterName("standard deviation",true);
  SourceTypeCmd->SetDefaultValue(1.);
  SourceTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

void MedLinacPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == EnergyCmd )
    { 
      G4double newEnergy = EnergyCmd->GetNewDoubleValue(newValue);
      MedLinacAction->SetEnergy(newEnergy);
 G4cout << "Energy "<< newEnergy << G4endl; 
    }

  if( command == SourceTypeCmd )
 {
   G4double newSourceType = SourceTypeCmd->GetNewDoubleValue(newValue);
      MedLinacAction->SetSourceType(newSourceType);
    
      G4cout << " standard deviation"<< newSourceType << G4endl;  
 }
}


MedLinacPrimaryGeneratorMessenger::~MedLinacPrimaryGeneratorMessenger()
{
  delete EnergyCmd;
  delete SourceTypeCmd;
}
