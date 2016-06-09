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
// Code developed by:
//  S.Guatelli
//
// $Id: BrachyRunMessenger.cc,v 1.4 2003/05/22 17:20:44 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// 
//
//    *******************************
//    *                             *
//    *    BrachyRunMessenger.cc    *
//    *                             *
//    *******************************

#include "BrachyRunMessenger.hh"
#include "BrachyRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

BrachyRunMessenger::BrachyRunMessenger( BrachyRunAction* pBrachyRun):
 runManager(pBrachyRun)
{ 
  runDir = new G4UIdirectory("/run/");
  runDir->SetGuidance("Control gamma energy.");
  
  primaryParticleEnergySpectrumCmd = 
                              new G4UIcmdWithAString("/run/energy",this);
  primaryParticleEnergySpectrumCmd->SetGuidance("Select the energy of gamma emitted by the source.");
  primaryParticleEnergySpectrumCmd->SetGuidance(" Iodium:  Iodium Source ");
  primaryParticleEnergySpectrumCmd->SetGuidance("  Iridium: Iridium  Source  "  );
  primaryParticleEnergySpectrumCmd->SetParameterName("choice",true);
  primaryParticleEnergySpectrumCmd->SetDefaultValue("Iridium");
  primaryParticleEnergySpectrumCmd->SetCandidates("Iridium / Iodium");
  primaryParticleEnergySpectrumCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

BrachyRunMessenger::~BrachyRunMessenger()
{
  delete primaryParticleEnergySpectrumCmd; 
  delete runDir;
}

void BrachyRunMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
 if(command == primaryParticleEnergySpectrumCmd)
  {
   if (newValue == "Iodium") runManager->SelectEnergy(1);
   else  runManager->SelectEnergy(0);
  }   
}

