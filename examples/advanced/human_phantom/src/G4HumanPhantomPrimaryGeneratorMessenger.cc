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
// 

#include "G4HumanPhantomPrimaryGeneratorMessenger.hh"
#include "G4HumanPhantomPrimaryGeneratorAction.hh"

#include "G4UIcmdWithAString.hh"


G4HumanPhantomPrimaryGeneratorMessenger::G4HumanPhantomPrimaryGeneratorMessenger(G4HumanPhantomPrimaryGeneratorAction* primaryGun)
  :primary(primaryGun)
{ 
  beamCmd = new G4UIcmdWithAString("/gun/setBeam",this);
  beamCmd->SetGuidance("Choose the type of beam");
  beamCmd->SetGuidance("Choice : beamAlongX, beamAlongY, beamAlongZ, isotropicFlux ");
  beamCmd->SetParameterName("choice",true);
  beamCmd->SetCandidates("beamAlongX beamAlongY beamAlongZ isotropicFlux");
  beamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

G4HumanPhantomPrimaryGeneratorMessenger::~G4HumanPhantomPrimaryGeneratorMessenger()
{
  delete beamCmd;
}

void G4HumanPhantomPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == beamCmd )
   { primary->SetBeam(newValue);
   }
}


