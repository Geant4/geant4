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


#include "ML2MainMessenger.hh"


CML2MainMessenger::CML2MainMessenger(CML2CInputData *CInputData)  
{
	this->CInputData=CInputData;
	this->phaseSpaceCentre=new G4UIcmdWith3VectorAndUnit("/general/centrePhaseSpace", this);
	this->phaseSpaceCentre->SetDefaultUnit("mm");
	this->phaseSpaceCentre->SetDefaultValue(G4ThreeVector(0.,0.,470.));

	this->phaseSpaceHalfSize=new G4UIcmdWith3VectorAndUnit("/general/halfSizePhaseSpace", this);
	this->phaseSpaceHalfSize->SetDefaultUnit("mm");
	this->phaseSpaceHalfSize->SetDefaultValue(G4ThreeVector(100.,100.,1.));

	this->bSavePhaseSpace=new G4UIcmdWithABool("/general/bSavePhaseSpace",this);
	this->bSavePhaseSpace->SetDefaultValue(false);

	this->bStopAtPhaseSpace=new G4UIcmdWithABool("/general/bStopAtPhaseSpace",this);
	this->bStopAtPhaseSpace->SetDefaultValue(false);

	this->bSaveROG=new G4UIcmdWithABool("/general/bSaveROG",this);
	this->bSaveROG->SetDefaultValue(false);

	this->bOnlyVisio=new G4UIcmdWithABool("/OnlyVisio",this);
	this->bOnlyVisio->SetDefaultValue(false);

	this->ROGOutFile=new G4UIcmdWithAString("/general/ROGOutFile",this);
	this->ROGOutFile->SetDefaultValue("");

	this->phaseSPaceOutFile=new G4UIcmdWithAString("/general/PhaseSpaceOutFile",this);
	this->phaseSPaceOutFile->SetDefaultValue("");


	this->minNumberOfEvents=new G4UIcmdWithAnInteger("/convergence/minNumberOfEvents", this);
	this->minNumberOfEvents->SetDefaultValue(10);

	this->bCompareExp=new G4UIcmdWithABool("/convergence/bCompareExp",this);
	this->bCompareExp->SetDefaultValue(false);

	this->fileExperimentalData=new G4UIcmdWithAString("/convergence/fileExperimentalData",this);
	this->fileExperimentalData->SetDefaultValue("");

	this->nBeam=new G4UIcmdWithAnInteger("/general/nBeam",this);
	this->nBeam->SetDefaultValue(0);
	
	this->nMaxParticlesInRamPlanePhaseSpace=new G4UIcmdWithAnInteger("/general/nMaxParticlesInRamPlanePhaseSpace",this);
	this->nMaxParticlesInRamPlanePhaseSpace->SetDefaultValue(0);

	this->saving_in_Selected_Voxels_every_events=new G4UIcmdWithAnInteger("/general/saving_in_Selected_Voxels_every_events",this);
	this->saving_in_Selected_Voxels_every_events->SetDefaultValue(1000);

	this->saving_in_ROG_Voxels_every_events=new G4UIcmdWithAnInteger("/general/saving_in_ROG_Voxels_every_events",this);
	this->saving_in_ROG_Voxels_every_events->SetDefaultValue(1000);

	this->max_N_particles_in_PhSp_File=new G4UIcmdWithAnInteger("/general/max_N_particles_in_PhSp_File",this);
	this->max_N_particles_in_PhSp_File->SetDefaultValue(1000);
}

CML2MainMessenger::~CML2MainMessenger(void)
{
	delete saving_in_Selected_Voxels_every_events; 
	delete saving_in_ROG_Voxels_every_events;
	delete max_N_particles_in_PhSp_File;
	delete phaseSpaceCentre; 
	delete phaseSpaceHalfSize;
	delete phaseSPaceOutFile;
	delete ROGOutFile;
	delete bSavePhaseSpace;
	delete bStopAtPhaseSpace;
	delete bSaveROG;

	delete nBeam;
	delete nMaxParticlesInRamPlanePhaseSpace;

	delete minNumberOfEvents;
	delete bCompareExp;
	delete bOnlyVisio;
	delete fileExperimentalData;
}

void CML2MainMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
	if (cmd==this->nBeam)
	{this->CInputData->setNBeams(this->nBeam->GetNewIntValue(newValue));}

	if (cmd==this->nMaxParticlesInRamPlanePhaseSpace)
	{this->CInputData->setNMaxParticlesInRamPlanePhaseSpace(this->nMaxParticlesInRamPlanePhaseSpace->GetNewIntValue(newValue));}

	if (cmd==this->saving_in_Selected_Voxels_every_events)
	{this->CInputData->setSaving_in_Selected_Voxels_every_events(this->saving_in_Selected_Voxels_every_events->GetNewIntValue(newValue));}

	if (cmd==this->saving_in_ROG_Voxels_every_events)
	{this->CInputData->setSaving_in_ROG_Voxels_every_events(this->saving_in_ROG_Voxels_every_events->GetNewIntValue(newValue));}

	if (cmd==this->max_N_particles_in_PhSp_File)
	{this->CInputData->setMax_N_particles_in_PhSp_File(this->max_N_particles_in_PhSp_File->GetNewIntValue(newValue));}

	if (cmd==this->phaseSpaceCentre)
	{this->CInputData->setPhaseSpaceCentre(this->phaseSpaceCentre->GetNew3VectorRawValue(newValue));}

	if (cmd==this->phaseSpaceHalfSize)
	{
		this->CInputData->setPhaseSpaceHalfSize(this->phaseSpaceHalfSize->GetNew3VectorRawValue(newValue));
	}

	if (cmd==this->bOnlyVisio)
	{this->CInputData->setbOnlyVisio(this->bOnlyVisio->GetNewBoolValue(newValue));}
	if (cmd==this->bSavePhaseSpace)
	{this->CInputData->setbSavePhaseSPace(this->bSavePhaseSpace->GetNewBoolValue(newValue));}
	if (cmd==this->bSaveROG)
	{this->CInputData->setbSaveROG(this->bSaveROG->GetNewBoolValue(newValue));}
	if (cmd==this->bStopAtPhaseSpace)
	{this->CInputData->setbStopAtPhaseSpace(this->bStopAtPhaseSpace->GetNewBoolValue(newValue));}
	if (cmd==this->phaseSPaceOutFile)
	{this->CInputData->setPhaseSpaceOutFile(newValue);}
	if (cmd==this->ROGOutFile)
	{this->CInputData->setROGOutFile(newValue);}
	if (cmd==this->minNumberOfEvents)
	{this->CInputData->setMinNumberOfEvents(this->minNumberOfEvents->GetNewIntValue(newValue));}
	if (cmd==this->bCompareExp)
	{this->CInputData->setBCompareExp(this->bCompareExp->GetNewBoolValue(newValue));}
	if (cmd==this->fileExperimentalData)
	{this->CInputData->setFileExperimentalData(newValue);}
}



