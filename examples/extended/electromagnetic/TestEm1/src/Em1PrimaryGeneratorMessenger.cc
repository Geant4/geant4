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
// $Id: Em1PrimaryGeneratorMessenger.cc,v 1.1 2001-12-07 11:49:10 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em1PrimaryGeneratorMessenger.hh"

#include "Em1PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em1PrimaryGeneratorMessenger::Em1PrimaryGeneratorMessenger(
                                             Em1PrimaryGeneratorAction* Em1Gun)
:Em1Action(Em1Gun)
{ 
  DefaultCmd = new G4UIcmdWithAnInteger("/gun/setDefault",this);
  DefaultCmd->SetGuidance("set/reset kinematic defined in PrimaryGenerator");
  DefaultCmd->SetGuidance("0=boxCenter, else=frontFace");
  DefaultCmd->SetParameterName("position",true);
  DefaultCmd->SetDefaultValue(0);
  DefaultCmd->AvailableForStates(PreInit,Idle);
  
  RndmCmd = new G4UIcmdWithADouble("/gun/rndm",this);
  RndmCmd->SetGuidance("random lateral extension on the beam");
  RndmCmd->SetGuidance("in fraction of 0.5*sizeYZ");
  RndmCmd->SetParameterName("rBeam",false);
  RndmCmd->SetRange("rBeam>=0.&&rBeam<=1.");
  RndmCmd->AvailableForStates(Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em1PrimaryGeneratorMessenger::~Em1PrimaryGeneratorMessenger()
{
  delete DefaultCmd;
  delete RndmCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{ 
  if (command == DefaultCmd)
   {Em1Action->SetDefaultKinematic(DefaultCmd->GetNewIntValue(newValue));}
   
  if (command == RndmCmd)
   {Em1Action->SetRndmBeam(RndmCmd->GetNewDoubleValue(newValue));}   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

