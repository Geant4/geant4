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


#include "ML2PhantomConstructionMessenger.hh"
#include "ML2PhantomConstruction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


CML2PhantomConstructionMessenger::CML2PhantomConstructionMessenger(CML2PhantomConstruction *phantomConstructor) : pPhantomConstructor (phantomConstructor)
{
	this->PhantomName=new G4UIcmdWithAString("/phantom/PhantomName",this);
	this->PhantomName->SetDefaultValue("myPhantom");

	this->PhantomSpecficationsFileName=new G4UIcmdWithAString("/phantom/PhantomSpecficationsFileName",this);
	this->PhantomSpecficationsFileName->SetDefaultValue("myPhantom");

	this->Phantom_nVoxelsX=new G4UIcmdWithAnInteger("/phantom/Phantom_nVoxelsX",this);
	this->Phantom_nVoxelsX->SetDefaultValue(100);

	this->Phantom_nVoxelsY=new G4UIcmdWithAnInteger("/phantom/Phantom_nVoxelsY",this);
	this->Phantom_nVoxelsY->SetDefaultValue(100);

	this->Phantom_nVoxelsZ=new G4UIcmdWithAnInteger("/phantom/Phantom_nVoxelsZ",this);
	this->Phantom_nVoxelsZ->SetDefaultValue(100);

	this->rotationX=new G4UIcmdWithADoubleAndUnit("/phantom/rotationX", this);
	this->rotationX->SetDefaultUnit("deg");
	this->rotationX->SetDefaultValue(0);

	this->rotationY=new G4UIcmdWithADoubleAndUnit("/phantom/rotationY", this);
	this->rotationY->SetDefaultUnit("deg");
	this->rotationY->SetDefaultValue(0);

	this->rotationZ=new G4UIcmdWithADoubleAndUnit("/phantom/rotationZ", this);
	this->rotationZ->SetDefaultUnit("deg");
	this->rotationZ->SetDefaultValue(0);
}

CML2PhantomConstructionMessenger::~CML2PhantomConstructionMessenger(void)
{
	delete PhantomName;
	delete Phantom_nVoxelsZ;
	delete Phantom_nVoxelsZ;
	delete Phantom_nVoxelsZ;
	delete rotationX;
	delete rotationY;
	delete rotationZ;
}
void CML2PhantomConstructionMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
	if (cmd==this->Phantom_nVoxelsX)
	{this->pPhantomConstructor->setPhantom_nVoxelsX(this->Phantom_nVoxelsX->GetNewIntValue(newValue));}

	if (cmd==this->Phantom_nVoxelsY)
	{this->pPhantomConstructor->setPhantom_nVoxelsY(this->Phantom_nVoxelsY->GetNewIntValue(newValue));}

	if (cmd==this->Phantom_nVoxelsZ)
	{this->pPhantomConstructor->setPhantom_nVoxelsZ(this->Phantom_nVoxelsZ->GetNewIntValue(newValue));}

	if (cmd==this->PhantomName)
	{this->pPhantomConstructor->setPhantomName(newValue);}

	if (cmd==this->PhantomSpecficationsFileName)
	{this->pPhantomConstructor->setPhantomSpecficationsFileName(newValue);}

	if (cmd==this->rotationX)
	{this->pPhantomConstructor->setPhantomRotationX(this->rotationX->GetNewDoubleValue(newValue));}

	if (cmd==this->rotationY)
	{this->pPhantomConstructor->setPhantomRotationY(this->rotationY->GetNewDoubleValue(newValue));}

	if (cmd==this->rotationZ)
	{this->pPhantomConstructor->setPhantomRotationZ(this->rotationZ->GetNewDoubleValue(newValue));}
}

