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
// $Id: Em5PrimaryGeneratorMessenger.cc,v 1.9 2002-12-16 16:30:08 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em5PrimaryGeneratorMessenger.hh"

#include "Em5PrimaryGeneratorAction.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em5PrimaryGeneratorMessenger::Em5PrimaryGeneratorMessenger(
                                             Em5PrimaryGeneratorAction* Em5Gun)
:Em5Action(Em5Gun)
{ 
 DefaultCmd = new G4UIcmdWithoutParameter("/testem/gun/setDefault",this);
 DefaultCmd->SetGuidance("set/reset the kinematic defined in PrimaryGenerator");
 DefaultCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  RndmCmd = new G4UIcmdWithADouble("/testem/gun/rndm",this);
  RndmCmd->SetGuidance("random lateral extension on the beam");
  RndmCmd->SetGuidance("in fraction of 0.5*sizeYZ");
  RndmCmd->SetParameterName("rBeam",false);
  RndmCmd->SetRange("rBeam>=0.&&rBeam<=1.");
  RndmCmd->AvailableForStates(G4State_Idle);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em5PrimaryGeneratorMessenger::~Em5PrimaryGeneratorMessenger()
{
  delete DefaultCmd;
  delete RndmCmd;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em5PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{ 
  if (command == DefaultCmd)
    {Em5Action->SetDefaultKinematic();}
   
  if (command == RndmCmd)
   { Em5Action->SetRndmBeam(RndmCmd->GetNewDoubleValue(newValue));}       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

