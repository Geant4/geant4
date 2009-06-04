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
// Code developed by:  S.Guatelli, guatelli@ge.infn.it
//
//    *****************************************
//    *                                       *
//    *   RadmonGunMessenger.cc *
//    *                                       *
//    *****************************************
//
//
// $Id: RadmonGunMessenger.cc,v 1.3 2009-06-04 13:45:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "RadmonGunMessenger.hh"
#include "RadmonPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

RadmonGunMessenger::RadmonGunMessenger(RadmonPrimaryGeneratorAction* prim): primary(prim)
{  
  gunDir = new G4UIdirectory("/gun/");
  gunDir->SetGuidance("Select the generation configuration of primary particles.");
 
  dataCmd = new G4UIcmdWithAString("/gun/data",this);
  dataCmd -> SetGuidance("Primary particle spectrum"); 
  dataCmd -> SetParameterName("choice",true);
  dataCmd -> AvailableForStates(G4State_PreInit,G4State_Idle); 

 //  uniformCmd = new G4UIcmdWithAString("/gun/uniform distribution",this);
//   uniformCmd -> SetGuidance("Generator momentum uniform"); 
//   uniformCmd -> SetParameterName("choice",true);
//   uniformCmd-> AvailableForStates(G4State_PreInit,G4State_Idle);

 }




RadmonGunMessenger::~RadmonGunMessenger()
{
  delete dataCmd;
  //  delete uniformCmd;
  delete gunDir;
} 
 
void RadmonGunMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  // if (command == dataCmd) primary -> ReadProbability(newValue);
 //if (command == dataCmd) primary -> ReadProbability(newValue);
  

}

