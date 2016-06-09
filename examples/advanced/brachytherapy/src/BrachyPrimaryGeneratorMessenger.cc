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
// $Id$
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
   // Change the energy of the emitted photons (iridium - iodium source)
  if(command == primaryParticleEnergySpectrumCmd)
    primaryAction -> SwitchEnergy(newValue); 
}

