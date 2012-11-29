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
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
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
	nRecycling=new G4UIcmdWithAnInteger("/primaryParticleData/nIdenticalParticles",this);
	nRecycling->SetDefaultValue(1);
	nRecycling->SetGuidance("number of identical particles generated in the primary generator");

	calculatedPhaseSpaceFileIN=new G4UIcmdWithAString("/primaryParticleData/calculatedPhaseSpaceFileIN",this);
	calculatedPhaseSpaceFileIN->SetDefaultValue("");
	calculatedPhaseSpaceFileIN->SetGuidance("full path and file name of the phase space file to be used as particle generator");

	sourceTypeName=new G4UIcmdWithAString("/primaryParticleData/sourceTypeName",this);
	sourceTypeName->SetDefaultValue("");
	sourceTypeName->SetGuidance("type of particle generator source  (randomTarget, phaseSpace)");

	nMaxParticlesInRamPhaseSpace=new G4UIcmdWithAnInteger("/primaryParticleData/nMaxParticlesInRamPhaseSpace",this);
	nMaxParticlesInRamPhaseSpace->SetDefaultValue(10000);
	nMaxParticlesInRamPhaseSpace->SetGuidance("maximum particle number loaded from the phase space file each time");

	GunMeanEnegy=new G4UIcmdWithADoubleAndUnit("/primaryParticleData/GunMeanEnegy", this);
	GunMeanEnegy->SetDefaultUnit("MeV");
	GunMeanEnegy->SetDefaultValue(6.);
	GunMeanEnegy->SetGuidance("mean energy of the primary particles");

	GunStdEnegy=new G4UIcmdWithADoubleAndUnit("/primaryParticleData/GunStdEnegy", this);
	GunStdEnegy->SetDefaultUnit("MeV");
	GunStdEnegy->SetDefaultValue(0.127);
	GunStdEnegy->SetGuidance("std energy of the primary particles");

	GunRadious=new G4UIcmdWithADoubleAndUnit("/primaryParticleData/GunRadious", this);
	GunRadious->SetDefaultUnit("mm");
	GunRadious->SetDefaultValue(10.);
	GunRadious->SetGuidance("radious primary particles beam");
}

CML2PrimaryGenerationActionMessenger::~CML2PrimaryGenerationActionMessenger(void)
{
	delete nRecycling;
	delete nMaxParticlesInRamPhaseSpace;
	delete GunMeanEnegy;
	delete GunStdEnegy;
	delete GunRadious;
	delete calculatedPhaseSpaceFileIN;
	delete sourceTypeName;
}
void CML2PrimaryGenerationActionMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
	if (cmd==GunMeanEnegy)
	{
		GunMeanEnegy->GetNewUnitValue(newValue);
		pML2PrimaryGenerationAction->setGunMeanEnergy(GunMeanEnegy->GetNewDoubleValue(newValue));
	}

	if (cmd==GunStdEnegy)
	{
		GunStdEnegy->GetNewUnitValue(newValue);
		pML2PrimaryGenerationAction->setGunStdEnergy(GunStdEnegy->GetNewDoubleValue(newValue));
	}

	if (cmd==GunRadious)
	{
		GunRadious->GetNewUnitValue(newValue);
		pML2PrimaryGenerationAction->setGunRadious(GunRadious->GetNewDoubleValue(newValue));
	}


	if (cmd==nMaxParticlesInRamPhaseSpace)
	{
		pML2PrimaryGenerationAction->setNMaxParticlesInRamPhaseSpace(nMaxParticlesInRamPhaseSpace->GetNewIntValue(newValue));
	}


	if (cmd==nRecycling)
	{pML2PrimaryGenerationAction->setNRecycling(nRecycling->GetNewIntValue(newValue));}

	if (cmd==calculatedPhaseSpaceFileIN)
	{pML2PrimaryGenerationAction->setCalculatedPhaseSpaceFileIN(newValue);}

	if (cmd==sourceTypeName)
	{pML2PrimaryGenerationAction->setSourceTypeName(newValue);}


}
