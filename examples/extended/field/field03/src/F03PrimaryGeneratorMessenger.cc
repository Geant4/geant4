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
// $Id: F03PrimaryGeneratorMessenger.cc,v 1.2 2001-07-11 09:58:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F03PrimaryGeneratorMessenger.hh"

#include "F03PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F03PrimaryGeneratorMessenger::F03PrimaryGeneratorMessenger(F03PrimaryGeneratorAction* F03Gun)
:F03Action(F03Gun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on, off(default)");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("off");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(PreInit,Idle);
 
  setxvertexCmd = new G4UIcmdWithADoubleAndUnit("/gun/xvertex",this);
  setxvertexCmd->SetGuidance(" Set x coord. of the primary vertex.");
  setxvertexCmd->SetParameterName("xv",true);
  setxvertexCmd->SetDefaultValue(0.0*mm) ; 
  
  setyvertexCmd = new G4UIcmdWithADoubleAndUnit("/gun/yvertex",this);
  setyvertexCmd->SetGuidance(" Set y coord. of the primary vertex.");
  setyvertexCmd->SetParameterName("yv",true);
  setyvertexCmd->SetDefaultValue(0.0*mm) ; 
  
  setzvertexCmd = new G4UIcmdWithADoubleAndUnit("/gun/zvertex",this);
  setzvertexCmd->SetGuidance(" Set z coord. of the primary vertex.");
  setzvertexCmd->SetParameterName("zv",true);
  setzvertexCmd->SetDefaultValue(0.0*mm) ; 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F03PrimaryGeneratorMessenger::~F03PrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete setxvertexCmd;
  delete setyvertexCmd;
  delete setzvertexCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F03PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmCmd )
   { F03Action->SetRndmFlag(newValue);}
  if( command == setxvertexCmd)
   { F03Action->Setxvertex(setxvertexCmd->GetNewDoubleValue(newValue));}
  if( command == setyvertexCmd)
   { F03Action->Setyvertex(setyvertexCmd->GetNewDoubleValue(newValue));}
  if( command == setzvertexCmd)
   { F03Action->Setzvertex(setzvertexCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

