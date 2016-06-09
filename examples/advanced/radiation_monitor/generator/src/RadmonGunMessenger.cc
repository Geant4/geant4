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
// Code developed by:  S.Guatelli, guatelli@ge.infn.it
//
//    *****************************************
//    *                                       *
//    *   RadmonGunMessenger.cc *
//    *                                       *
//    *****************************************
//
//
// $Id: RadmonGunMessenger.cc,v 1.3.2.2 2009/08/11 14:20:35 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-04 $
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

