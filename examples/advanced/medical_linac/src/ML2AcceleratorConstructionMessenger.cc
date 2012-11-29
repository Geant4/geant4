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

#include "ML2AcceleratorConstructionMessenger.hh"
#include "ML2AcceleratorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

CML2AcceleratorConstructionMessenger::CML2AcceleratorConstructionMessenger(CML2AcceleratorConstruction *acceleratorConstructor) : pAcceleratorConstructor (acceleratorConstructor)
{
	bOnlyVisio=false;
	AcceleratorName=new G4UIcmdWithAString("/accelerator/AcceleratorName",this);
	AcceleratorName->SetDefaultValue("acc1");
	AcceleratorName->SetGuidance("accelerator name to select among those implemented (acc1)");
	pAcceleratorConstructor->setAcceleratorName("acc1");

	acceleratorMacFileName=new G4UIcmdWithAString("/accelerator/AcceleratorMacFileName",this);
	acceleratorMacFileName->SetDefaultValue("");
	acceleratorMacFileName->SetGuidance("full path and macro file name containing specific setup data for the accelerator chosen");
	pAcceleratorConstructor->setAcceleratorMacFileName("");


	rotationX  =new G4UIcmdWithADoubleAndUnit("/accelerator/rotationX", this);
	rotationX  ->SetDefaultUnit("deg");
	rotationX  ->SetDefaultValue(0.);
	rotationX->SetGuidance("angles of rotation along X [deg]");

	isoCentre=new G4UIcmdWithADoubleAndUnit("/accelerator/isoCentre", this);
	isoCentre->SetDefaultUnit("mm");
	isoCentre->SetDefaultValue(1000.);
	isoCentre->SetGuidance("distance between the isocentre and the target of the accelerator");
	pAcceleratorConstructor->setIsoCentre(1000.*mm);


	bRotate90Y  =new G4UIcmdWithABool("/accelerator/rotation90Y", this);
	bRotate90Y  ->SetDefaultValue(false);
	bRotate90Y->SetGuidance("to rotate the accelerator of 90 deg around the Y axis (true)");
	pAcceleratorConstructor->setRotation90Y(false);
}

CML2AcceleratorConstructionMessenger::~CML2AcceleratorConstructionMessenger(void)
{
	delete AcceleratorName;
	delete rotationX ;
	delete acceleratorMacFileName;
	delete isoCentre;
	delete bRotate90Y;
}
void CML2AcceleratorConstructionMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{

	if (cmd==AcceleratorName)
	{pAcceleratorConstructor->setAcceleratorName(newValue);}

	if (cmd==acceleratorMacFileName)
	{pAcceleratorConstructor->setAcceleratorMacFileName(newValue);}


	if (cmd==rotationX)
	{
		if (bOnlyVisio)
		{
			G4RotationMatrix *rm=pAcceleratorConstructor->rotateAccelerator(rotationX ->GetNewDoubleValue(newValue));
			CML2PrimaryGenerationAction::GetInstance()->setRotation(rm);
			CML2PhantomConstruction::GetInstance()->resetSensDet();
// what follows seems to be necessary to have a good refresh
			G4UImanager* UI = G4UImanager::GetUIpointer();
			G4String command;
			command = "/run/beamOn 0";
			UI->ApplyCommand(command); 
			command = "/vis/viewer/flush";
			UI->ApplyCommand(command); 
		}
		else
		{
			pAcceleratorConstructor->addAcceleratorRotationsX(rotationX ->GetNewDoubleValue(newValue));
		}
	}

	if (cmd==isoCentre)
	{
		isoCentre->GetNewUnitValue(newValue);
		pAcceleratorConstructor->setIsoCentre(isoCentre->GetNewDoubleValue(newValue));
	}
	if (cmd==bRotate90Y)
	{pAcceleratorConstructor->setRotation90Y(bRotate90Y->GetNewBoolValue(newValue));}
}
