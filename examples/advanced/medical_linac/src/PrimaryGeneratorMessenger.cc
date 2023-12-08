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
// Code developed by
// Silvia Pozzi (1), silvia.pozzi@iss.it
// Barbara Caccia (1), barbara.caccia@iss.it
// Carlo Mancini Terracciano (2), carlo.mancini.terracciano@roma1.infn.it
// (1) Istituto Superiore di Sanita' and INFN Roma, Italy
// (2) Univ. La Sapienza and INFN Roma, Italy

#include "PrimaryGeneratorMessenger.hh"

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete cmdGunRadius;
  delete cmdGunMeanEnergy;
  delete cmdGunStdEnergy;
}

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* myGun)
: primaryGenPointer(myGun), cmdGunRadius(nullptr),
  cmdGunMeanEnergy(nullptr), cmdGunStdEnergy(nullptr)
{
  cmdGunRadius = new G4UIcmdWithADoubleAndUnit("/PrimaryGenerator/gunRadius", this);
  cmdGunRadius -> SetDefaultUnit("mm");
  cmdGunRadius -> SetDefaultValue(1.);
  primaryGenPointer -> SetGunRadius(1.*mm);
  
  cmdGunMeanEnergy = new G4UIcmdWithADoubleAndUnit("/PrimaryGenerator/gunMeanEnergy", this);
  cmdGunMeanEnergy -> SetDefaultUnit("MeV");
  cmdGunMeanEnergy -> SetDefaultValue(12.);
  primaryGenPointer -> SetGunMeanEnergy(12.*MeV);
  
  cmdGunStdEnergy = new G4UIcmdWithADoubleAndUnit("/PrimaryGenerator/gunStdEnergy", this);
  cmdGunStdEnergy -> SetDefaultUnit("MeV");
  cmdGunStdEnergy -> SetDefaultValue(0.127);
  primaryGenPointer -> SetGunStdEnergy(0.127*MeV);
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
  if ( cmd == cmdGunRadius )
	{
		cmdGunRadius -> GetNewUnitValue(newValue);
		primaryGenPointer -> SetGunRadius(cmdGunRadius -> GetNewDoubleValue(newValue));
	}
  else if ( cmd == cmdGunMeanEnergy )
	{
		cmdGunMeanEnergy -> GetNewUnitValue(newValue);
		primaryGenPointer -> SetGunMeanEnergy(cmdGunMeanEnergy -> GetNewDoubleValue(newValue));
	}
  else if ( cmd == cmdGunStdEnergy )
	{
		cmdGunStdEnergy -> GetNewUnitValue(newValue);
		primaryGenPointer -> SetGunStdEnergy(cmdGunStdEnergy -> GetNewDoubleValue(newValue));
	}
  else
    G4cerr << "PrimaryGeneratorMessenger::SetNewValue: command not found" << G4endl;
  
}
