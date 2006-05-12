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
//    *****************************************
//    *                                       *
//    *    BrachyPrimaryGeneratorMessenger.cc *
//    *                                       *
//    *****************************************
//
//
// $Id: BrachyPrimaryGeneratorMessenger.cc,v 1.1 2006-05-12 17:08:06 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "BrachyPrimaryGeneratorMessenger.hh"
#include "BrachyPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

BrachyPrimaryGeneratorMessenger::BrachyPrimaryGeneratorMessenger(BrachyPrimaryGeneratorAction* primary):
 primaryAction(primary)
{ 
  primaryDir = new G4UIdirectory("/primary/");
  primaryDir -> SetGuidance("primary particles control.");
 
  // Define interactive command
  primaryParticleEnergySpectrumCmd = 
                              new G4UIcmdWithAString("/primary/energy",this);
  primaryParticleEnergySpectrumCmd -> SetGuidance("Select the energy of gamma emitted by the source.");
  primaryParticleEnergySpectrumCmd -> SetGuidance("Iodium: Iodium Source");
  primaryParticleEnergySpectrumCmd -> SetGuidance("Iridium: Iridium Source");
  primaryParticleEnergySpectrumCmd -> SetParameterName("choice",true);
  primaryParticleEnergySpectrumCmd -> SetDefaultValue("Iridium");
  primaryParticleEnergySpectrumCmd -> SetCandidates("Iridium / Iodium");
  primaryParticleEnergySpectrumCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);      
}

BrachyPrimaryGeneratorMessenger::~BrachyPrimaryGeneratorMessenger()
{
  delete primaryParticleEnergySpectrumCmd;
  delete primaryDir;
}

void BrachyPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
   // Change the energy of the emitted photons
  if(command == primaryParticleEnergySpectrumCmd)
    primaryAction -> SwitchEnergy(newValue); 
}

