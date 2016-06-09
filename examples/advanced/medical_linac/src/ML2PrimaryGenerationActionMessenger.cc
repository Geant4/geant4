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
// The code was written by :
//	^Claudio Andenna claudio.andenna@iss.infn.it, claudio.andenna@ispesl.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//
// ^ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2PrimaryGenerationActionMessenger.hh"
#include "ML2PrimaryGenerationAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

CML2PrimaryGenerationActionMessenger::CML2PrimaryGenerationActionMessenger(CML2PrimaryGenerationAction *PML2PrimaryGenerationAction) : pML2PrimaryGenerationAction(PML2PrimaryGenerationAction)
{
	this->nIdenticalParticles=new G4UIcmdWithAnInteger("/primaryParticleData/nIdenticalParticles",this);
	this->nIdenticalParticles->SetDefaultValue(1);

	this->calculatedPhaseSpaceFileIN=new G4UIcmdWithAString("/primaryParticleData/calculatedPhaseSpaceFileIN",this);
	this->calculatedPhaseSpaceFileIN->SetDefaultValue("");

	this->sourceTypeName=new G4UIcmdWithAString("/primaryParticleData/sourceTypeName",this);
	this->sourceTypeName->SetDefaultValue("");

	this->nMaxParticlesInRamPhaseSpace=new G4UIcmdWithAnInteger("/primaryParticleData/nMaxParticlesInRamPhaseSpace",this);
	this->nMaxParticlesInRamPhaseSpace->SetDefaultValue(10000);

	this->nLoopsPhSpParticles=new G4UIcmdWithAnInteger("/primaryParticleData/nLoopsPhSpParticles",this);
	this->nLoopsPhSpParticles->SetDefaultValue(1);

	this->GunMeanEnegy=new G4UIcmdWithADoubleAndUnit("/primaryParticleData/GunMeanEnegy", this);
	this->GunMeanEnegy->SetDefaultUnit("MeV");
	this->GunMeanEnegy->SetDefaultValue(6.);

	this->GunStdEnegy=new G4UIcmdWithADoubleAndUnit("/primaryParticleData/GunStdEnegy", this);
	this->GunStdEnegy->SetDefaultUnit("MeV");
	this->GunStdEnegy->SetDefaultValue(0.127);

	this->GunRadious=new G4UIcmdWithADoubleAndUnit("/primaryParticleData/GunRadious", this);
	this->GunRadious->SetDefaultUnit("mm");
	this->GunRadious->SetDefaultValue(10.);
}

CML2PrimaryGenerationActionMessenger::~CML2PrimaryGenerationActionMessenger(void)
{
	delete nIdenticalParticles;
	delete nLoopsPhSpParticles;
	delete nMaxParticlesInRamPhaseSpace;
	delete GunMeanEnegy;
	delete GunStdEnegy;
	delete GunRadious;
	delete calculatedPhaseSpaceFileIN;
	delete sourceTypeName;
}
void CML2PrimaryGenerationActionMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
	if (cmd==this->GunMeanEnegy)
	{
		this->GunMeanEnegy->GetNewUnitValue(newValue);
		this->pML2PrimaryGenerationAction->setGunMeanEnergy(this->GunMeanEnegy->GetNewDoubleValue(newValue));
	}

	if (cmd==this->GunStdEnegy)
	{
		this->GunStdEnegy->GetNewUnitValue(newValue);
		this->pML2PrimaryGenerationAction->setGunStdEnergy(this->GunStdEnegy->GetNewDoubleValue(newValue));
	}

	if (cmd==this->GunRadious)
	{
		this->GunRadious->GetNewUnitValue(newValue);
		this->pML2PrimaryGenerationAction->setGunRadious(this->GunRadious->GetNewDoubleValue(newValue));
	}


	if (cmd==this->nMaxParticlesInRamPhaseSpace)
	{
		this->pML2PrimaryGenerationAction->setNMaxParticlesInRamPhaseSpace(this->nMaxParticlesInRamPhaseSpace->GetNewIntValue(newValue));
	}


	if (cmd==this->nIdenticalParticles)
	{this->pML2PrimaryGenerationAction->setNIdenticalParticles(this->nIdenticalParticles->GetNewIntValue(newValue));}

	if (cmd==this->calculatedPhaseSpaceFileIN)
	{this->pML2PrimaryGenerationAction->setCalculatedPhaseSpaceFileIN(newValue);}

	if (cmd==this->sourceTypeName)
	{this->pML2PrimaryGenerationAction->setSourceTypeName(newValue);}


}
