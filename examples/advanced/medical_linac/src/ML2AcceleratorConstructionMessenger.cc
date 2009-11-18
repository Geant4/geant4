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


#include "ML2AcceleratorConstructionMessenger.hh"
#include "ML2AcceleratorConstruction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

CML2AcceleratorConstructionMessenger::CML2AcceleratorConstructionMessenger(CML2AcceleratorConstruction *acceleratorConstructor) : pAcceleratorConstructor (acceleratorConstructor)
{
	this->AcceleratorName=new G4UIcmdWithAString("/accelerator/AcceleratorName",this);
	this->AcceleratorName->SetDefaultValue("nome del mio acceleratore");

	this->acceleratorSpecficationsFile=new G4UIcmdWithAString("/accelerator/AcceleratorSpecficationsFileName",this);
	this->acceleratorSpecficationsFile->SetDefaultValue("");

	this->rotationX=new G4UIcmdWithADoubleAndUnit("/accelerator/rotationX", this);
	this->rotationX->SetDefaultUnit("deg");
	this->rotationX->SetDefaultValue(0);

	this->rotationY=new G4UIcmdWithADoubleAndUnit("/accelerator/rotationY", this);
	this->rotationY->SetDefaultUnit("deg");
	this->rotationY->SetDefaultValue(0);

	this->rotationZ=new G4UIcmdWithADoubleAndUnit("/accelerator/rotationZ", this);
	this->rotationZ->SetDefaultUnit("deg");
	this->rotationZ->SetDefaultValue(0);
}

CML2AcceleratorConstructionMessenger::~CML2AcceleratorConstructionMessenger(void)
{
	delete AcceleratorName;
	delete rotationX;
	delete rotationY;
	delete rotationZ;
	delete acceleratorSpecficationsFile;
}
void CML2AcceleratorConstructionMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{

	if (cmd==this->AcceleratorName)
	{this->pAcceleratorConstructor->setAcceleratorName(newValue);}

	if (cmd==this->acceleratorSpecficationsFile)
	{this->pAcceleratorConstructor->setAcceleratorSpecficationsFileName(newValue);}

	if (cmd==this->rotationX)
	{this->pAcceleratorConstructor->setAcceleratorRotationX(this->rotationX->GetNewDoubleValue(newValue));}

	if (cmd==this->rotationY)
	{this->pAcceleratorConstructor->setAcceleratorRotationY(this->rotationY->GetNewDoubleValue(newValue));}

	if (cmd==this->rotationZ)
	{this->pAcceleratorConstructor->setAcceleratorRotationZ(this->rotationZ->GetNewDoubleValue(newValue));}
}
