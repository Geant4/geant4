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
// Code developed by:
//  S.Guatelli
//
//    *********************************
//    *                               *
//    *    RemSimPrimaryGeneratorMessenger.cc *
//    *                               *
//    *********************************
//
//
// $Id: RemSimPrimaryGeneratorMessenger.cc,v 1.2 2004-03-12 10:55:55 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "RemSimPrimaryGeneratorMessenger.hh"
#include "RemSimRunAction.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

RemSimPrimaryGeneratorMessenger::RemSimPrimaryGeneratorMessenger( RemSimPrimaryGeneratorAction* prim): primary(prim)
{ 
  primaryDir = new G4UIdirectory("/beam/");
  primaryDir->SetGuidance("primary particle generator control.");
        
  isotropicFluxCmd = new G4UIcmdWithAString("/beam/flux",this);
  isotropicFluxCmd->SetGuidance("Assign the isotropic flux option to primary particles."); 
  isotropicFluxCmd->SetParameterName("choice",true);
  isotropicFluxCmd->SetCandidates("isotropic");
  isotropicFluxCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
 }

RemSimPrimaryGeneratorMessenger::~RemSimPrimaryGeneratorMessenger()
{
  delete isotropicFluxCmd;
  delete primaryDir;
}

void RemSimPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
   if(command == isotropicFluxCmd )
   {
     if(newValue=="isotropic") 
       primary->GenerateIsotropicFlux();
     else G4cout<< newValue<< " option is not available!"<<G4cout; 
  }
}

