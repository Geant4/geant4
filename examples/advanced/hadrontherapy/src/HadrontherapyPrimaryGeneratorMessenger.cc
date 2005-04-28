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
// $Id: HadrontherapyPrimaryGeneratorMessenger.cc,v 1.2 2005-04-28 20:39:33 mpiergen Exp $
//
//
// Code developed by: M. Piergentili
//
#include "HadrontherapyPrimaryGeneratorMessenger.hh"

#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

HadrontherapyPrimaryGeneratorMessenger::HadrontherapyPrimaryGeneratorMessenger(HadrontherapyPrimaryGeneratorAction* HadrontherapyGun)
  :HadrontherapyAction(HadrontherapyGun)

{  
  beamDir = new G4UIdirectory("/beam/");

  EnergyCmd = new G4UIcmdWithADoubleAndUnit("/beam/energy",this);
  EnergyCmd->SetGuidance("Select the energy (and unit) of the primary particles");
  EnergyCmd->SetGuidance("Type particle Energy");
  EnergyCmd->SetParameterName("energy",true,true);
  EnergyCmd->SetDefaultValue(1.);
  EnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SourceTypeCmd = new G4UIcmdWithADoubleAndUnit("/beam/energySpread",this);
  SourceTypeCmd->SetGuidance("Select the standard deviation of the energy.");
  SourceTypeCmd->SetParameterName("sigma",true);
  SourceTypeCmd->SetDefaultValue(1.);
  SourceTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BeamCmd = new G4UIcmdWithADoubleAndUnit("/beam/spotXPosition",this);
  BeamCmd->SetGuidance("Select the x position (and unit) of the spot");
  BeamCmd->SetParameterName("x0",true,true);
  BeamCmd->SetDefaultValue(1.);
  BeamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SizeYCmd = new G4UIcmdWithADoubleAndUnit("/beam/spotSizeY",this);
  SizeYCmd->SetGuidance("Select the size (and unit) of the spot");
  SizeYCmd->SetParameterName("sizeY",true,true);
  SizeYCmd->SetDefaultValue(1.);
  SizeYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SizeZCmd = new G4UIcmdWithADoubleAndUnit("/beam/spotSizeZ",this);
  SizeZCmd->SetGuidance("Select the size (and unit) of the spot");
  SizeZCmd->SetParameterName("sizeZ",true,true);
  SizeZCmd->SetDefaultValue(1.);
  SizeZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

void HadrontherapyPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == EnergyCmd )
    { 
      G4double newEnergy = EnergyCmd->GetNewDoubleValue(newValue);
      HadrontherapyAction->SetEnergy(newEnergy);
 G4cout << "Energy "<< newEnergy << G4endl; 
    }

  if( command == SourceTypeCmd )
 {
   G4double newSourceType = SourceTypeCmd->GetNewDoubleValue(newValue);
      HadrontherapyAction->SetSourceType(newSourceType);
    
      G4cout << " standard deviation"<< newSourceType << G4endl;  
 }

  if( command == BeamCmd )
 {
   G4double newPosition = BeamCmd->GetNewDoubleValue(newValue);
      HadrontherapyAction->SetXPosition(newPosition);
    
      G4cout << " beam source x position"<< newPosition << G4endl;  
 }

  if( command == SizeYCmd )
 {
   G4double newSpotSizeY = SizeYCmd->GetNewDoubleValue(newValue);
      HadrontherapyAction->SetSpotSizeY(newSpotSizeY);
    
      G4cout << " beam source y size"<< newSpotSizeY << G4endl;  
 }
  if( command == SizeZCmd )
 {
   G4double newSpotSizeZ = SizeZCmd->GetNewDoubleValue(newValue);
      HadrontherapyAction->SetSpotSizeZ(newSpotSizeZ);
    
      G4cout << " beam source y size"<< newSpotSizeZ << G4endl;  
 }

}


HadrontherapyPrimaryGeneratorMessenger::~HadrontherapyPrimaryGeneratorMessenger()
{
  delete EnergyCmd;
  delete SourceTypeCmd;
  delete beamDir;
  delete SizeYCmd;
  delete SizeZCmd;
}
