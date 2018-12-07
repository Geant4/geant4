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


#include "ML2MainMessenger.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

CML2MainMessenger::CML2MainMessenger(CML2CInputData *InData)  
{
	CInputData = InData;
	phaseSpaceCentre = new G4UIcmdWith3VectorAndUnit("/general/centrePhaseSpace", this);
	phaseSpaceCentre -> SetDefaultUnit("mm");
	phaseSpaceCentre -> SetDefaultValue(G4ThreeVector(0.,0.,164.));
	phaseSpaceCentre -> SetGuidance("position of the centre of the phase space plane");
	CInputData -> setPhaseSpaceCentre(G4ThreeVector(0.*mm,0.*mm,164.*mm));

	phaseSpaceHalfSize = new G4UIcmdWith3VectorAndUnit("/general/halfSizePhaseSpace", this);
	phaseSpaceHalfSize -> SetDefaultUnit("mm");
	phaseSpaceHalfSize -> SetDefaultValue(G4ThreeVector(100.,100.,1.));
	phaseSpaceHalfSize -> SetGuidance("half size of the phase space plane");
	CInputData -> setPhaseSpaceHalfSize(G4ThreeVector(100.*mm,100.*mm,1.*mm));

	bSavePhaseSpace = new G4UIcmdWithABool("/general/bSavePhaseSpace",this);
	bSavePhaseSpace -> SetDefaultValue(false);
	bSavePhaseSpace -> SetGuidance("save the phase space");
	CInputData -> setbSavePhaseSPace(false);

	bForcePhaseSpaceBeforeJaws = new G4UIcmdWithABool("/general/bForcePhaseSpaceBeforeJaws",this);
	bForcePhaseSpaceBeforeJaws -> SetDefaultValue(false);
	bForcePhaseSpaceBeforeJaws -> SetGuidance("automatically put the phase plane before the jaws");
	CInputData -> setbForcePhaseSpaceBeforeJaws(false);

	bStopAtPhaseSpace = new G4UIcmdWithABool("/general/bStopAtPhaseSpace",this);
	bStopAtPhaseSpace -> SetDefaultValue(false);
	bStopAtPhaseSpace -> SetGuidance("kill the particle at the phase space");
	CInputData-> setbStopAtPhaseSpace(false);

	bSaveROG = new G4UIcmdWithABool("/general/bSaveROG",this);
	bSaveROG -> SetDefaultValue(true);
	bSaveROG -> SetGuidance("save the ROG volume");
	CInputData -> setbSaveROG(true);

	bOnlyVisio = new G4UIcmdWithABool("/OnlyVisio",this);
	bOnlyVisio -> SetDefaultValue(false);
	bOnlyVisio -> SetGuidance("switch the visualization mode");
	CInputData -> setbOnlyVisio(false);

	ROGOutFile = new G4UIcmdWithAString("/general/ROGOutFile",this);
	ROGOutFile -> SetDefaultValue("");
	ROGOutFile -> SetGuidance("full path and name of the ROG file");
	CInputData -> setROGOutFile("defaultROGFile.txt");

	phaseSPaceOutFile = new G4UIcmdWithAString("/general/PhaseSpaceOutFile",this);
	phaseSPaceOutFile -> SetDefaultValue("");
	phaseSPaceOutFile -> SetGuidance("full path and name of the phase space file");
	CInputData -> setPhaseSpaceOutFile("");

	maxNumberOfEvents=new G4UIcmdWithAnInteger("/convergence/maxNumberOfEvents", this);
	maxNumberOfEvents->SetDefaultValue(10);
	maxNumberOfEvents->SetGuidance("number of events to be reached in one experimental voxel");
	CInputData->setMaxNumberOfEvents(10);

	nMaxLoop = new G4UIcmdWithAnInteger("/convergence/nMaxLoop", this);
	nMaxLoop -> SetDefaultValue(1);
	nMaxLoop -> SetGuidance("used if /convergence/bCompareExp is false");
	CInputData -> setNmaxLoop(1);

	bCompareExp = new G4UIcmdWithABool("/convergence/bCompareExp", this);
	bCompareExp -> SetDefaultValue(false);
	bCompareExp -> SetGuidance("compare the data with an experimental data file");
	CInputData -> setBCompareExp(false);

	fileExperimentalData = new G4UIcmdWithAString("/convergence/fileExperimentalData", this);
	fileExperimentalData -> SetDefaultValue("");
	fileExperimentalData -> SetGuidance("full path and name of the experimental data file");
	CInputData -> setFileExperimentalData("");

	fileExperimentalDataOut = new G4UIcmdWithAString("/convergence/fileExperimentalDataOut", this);
	fileExperimentalDataOut -> SetDefaultValue("");
	fileExperimentalDataOut -> SetGuidance("full path and name of the experimental data out file");
	CInputData -> setFileExperimentalDataOut("");

	nBeam = new G4UIcmdWithAnInteger("/general/nBeam", this);
	nBeam -> SetDefaultValue(100);
	nBeam -> SetGuidance("number of events to run");
	CInputData -> setNBeams(100);

	nMaxParticlesInRamPlanePhaseSpace = new G4UIcmdWithAnInteger("/general/nMaxParticlesInRamPlanePhaseSpace",this);
	nMaxParticlesInRamPlanePhaseSpace -> SetDefaultValue(10000);
	nMaxParticlesInRamPlanePhaseSpace -> SetGuidance("maximum number of particles stored in RAM before saving phase space file");
	CInputData -> setNMaxParticlesInRamPlanePhaseSpace(10000);

	saving_in_Selected_Voxels_every_events = new G4UIcmdWithAnInteger("/general/saving_in_Selected_Voxels_every_events",this);
	saving_in_Selected_Voxels_every_events -> SetDefaultValue(10000);
	saving_in_Selected_Voxels_every_events -> SetGuidance("maximum number of particles stored before saving experimental data out file");
	CInputData -> setSaving_in_Selected_Voxels_every_events(10000);

	saving_in_ROG_Voxels_every_events = new G4UIcmdWithAnInteger("/general/saving_in_ROG_Voxels_every_events",this);
	saving_in_ROG_Voxels_every_events -> SetDefaultValue(1000);
	saving_in_ROG_Voxels_every_events -> SetGuidance("maximum number of particles stored before saving ROG file");
	CInputData -> setSaving_in_ROG_Voxels_every_events(1000);

	max_N_particles_in_PhSp_File = new G4UIcmdWithAnInteger("/general/max_N_particles_in_PhSp_File",this);
	max_N_particles_in_PhSp_File -> SetDefaultValue(1000);
	max_N_particles_in_PhSp_File -> SetGuidance("maximum number of particles stored in the phase space file");
	CInputData -> setMax_N_particles_in_PhSp_File(1000);

       // SUSANNA: command to fix the number of voxels in the phantom segmentation
        fVoxelsXCmd = new G4UIcmdWithAnInteger("/general/numberVoxelsX",this);
	fVoxelsXCmd -> SetDefaultValue(30);
	fVoxelsXCmd -> SetGuidance("number of voxels along X axis");
	CInputData->setVoxelsX(30);

        fVoxelsYCmd = new G4UIcmdWithAnInteger("/general/numberVoxelsY",this);
	fVoxelsYCmd -> SetDefaultValue(30);
	fVoxelsYCmd -> SetGuidance("number of voxels along Y axis");
        CInputData->setVoxelsY(30);

        fVoxelsZCmd = new G4UIcmdWithAnInteger("/general/numberVoxelsZ",this);
	fVoxelsZCmd -> SetDefaultValue(30);
	fVoxelsZCmd -> SetGuidance("number of voxels along Z axis");
	CInputData->setVoxelsZ(30);
       
}

CML2MainMessenger::~CML2MainMessenger(void)
{
        delete fVoxelsXCmd;
        delete fVoxelsYCmd;
        delete fVoxelsZCmd;
	delete saving_in_Selected_Voxels_every_events; 
	delete saving_in_ROG_Voxels_every_events;
	delete max_N_particles_in_PhSp_File;
	delete phaseSpaceCentre; 
	delete phaseSpaceHalfSize;
	delete phaseSPaceOutFile;
	delete ROGOutFile;
	delete bSavePhaseSpace;
	delete bForcePhaseSpaceBeforeJaws;
	delete bStopAtPhaseSpace;
	delete bSaveROG;

	delete nBeam;
	delete nMaxParticlesInRamPlanePhaseSpace;

	delete maxNumberOfEvents;
	delete bCompareExp;
	delete bOnlyVisio;
	delete fileExperimentalData;
	delete fileExperimentalDataOut;
	delete nMaxLoop;
}

void CML2MainMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
	if (cmd==nBeam)
	{CInputData->setNBeams(nBeam->GetNewIntValue(newValue));}

	if (cmd==nMaxParticlesInRamPlanePhaseSpace)
	{CInputData->setNMaxParticlesInRamPlanePhaseSpace(nMaxParticlesInRamPlanePhaseSpace->GetNewIntValue(newValue));}

	if (cmd==saving_in_Selected_Voxels_every_events)
	{CInputData->setSaving_in_Selected_Voxels_every_events(saving_in_Selected_Voxels_every_events->GetNewIntValue(newValue));}

	if (cmd==saving_in_ROG_Voxels_every_events)
	{CInputData->setSaving_in_ROG_Voxels_every_events(saving_in_ROG_Voxels_every_events->GetNewIntValue(newValue));}

	if (cmd==max_N_particles_in_PhSp_File)
	{CInputData->setMax_N_particles_in_PhSp_File(max_N_particles_in_PhSp_File->GetNewIntValue(newValue));}

	if (cmd==phaseSpaceCentre)
	{CInputData->setPhaseSpaceCentre(phaseSpaceCentre->GetNew3VectorRawValue(newValue));}

	if (cmd==phaseSpaceHalfSize)
	{
		CInputData->setPhaseSpaceHalfSize(phaseSpaceHalfSize->GetNew3VectorRawValue(newValue));
	}

	if (cmd==bOnlyVisio)
	{CInputData->setbOnlyVisio(bOnlyVisio->GetNewBoolValue(newValue));}

	if (cmd==bForcePhaseSpaceBeforeJaws)
	{CInputData->setbForcePhaseSpaceBeforeJaws(bForcePhaseSpaceBeforeJaws->GetNewBoolValue(newValue));}

	if (cmd==bSavePhaseSpace)
	{CInputData->setbSavePhaseSPace(bSavePhaseSpace->GetNewBoolValue(newValue));}

	if (cmd==bSaveROG)
	{CInputData->setbSaveROG(bSaveROG->GetNewBoolValue(newValue));}

	if (cmd==bStopAtPhaseSpace)
	{CInputData->setbStopAtPhaseSpace(bStopAtPhaseSpace->GetNewBoolValue(newValue));}

	if (cmd==phaseSPaceOutFile)
	{CInputData->setPhaseSpaceOutFile(newValue);}

	if (cmd==ROGOutFile)
	{CInputData->setROGOutFile(newValue);}

	if (cmd==maxNumberOfEvents)
	{CInputData->setMaxNumberOfEvents(maxNumberOfEvents->GetNewIntValue(newValue));}

	if (cmd==nMaxLoop)
	{CInputData->setNmaxLoop(nMaxLoop->GetNewIntValue(newValue));}

	if (cmd==bCompareExp)
	{CInputData->setBCompareExp(bCompareExp->GetNewBoolValue(newValue));}

	if (cmd==fileExperimentalData)
	{CInputData->setFileExperimentalData(newValue);}

	if (cmd==fileExperimentalDataOut)
	{CInputData->setFileExperimentalDataOut(newValue);} 


        if (cmd==fVoxelsXCmd) {CInputData->setVoxelsX(fVoxelsXCmd->GetNewIntValue(newValue));}
        if (cmd==fVoxelsYCmd) {CInputData->setVoxelsY(fVoxelsYCmd->GetNewIntValue(newValue));}
       if (cmd==fVoxelsZCmd) {CInputData->setVoxelsZ(fVoxelsZCmd->GetNewIntValue(newValue));}
}
