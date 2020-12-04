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


#include "ML2Acc1Messenger.hh"
#include "ML2Acc1.hh"
#include "ML2Accelerator.hh"

#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

CML2Acc1Messenger::CML2Acc1Messenger(CML2Acc1 *acc1) : pAcc1(acc1)
{
	idEnergy=new G4UIcmdWithAnInteger("/accelerator/idEnergy",this);
	idEnergy->SetDefaultValue(6);
	pAcc1->setidEnergy(6);

	leavesA=new G4UIcmdWithADoubleAndUnit("/accelerator/leavesA", this);
	leavesA->SetDefaultUnit("mm");
	leavesA->SetDefaultValue(300.);
	pAcc1->setLeavesAx(300.*mm);

	leavesB=new G4UIcmdWithADoubleAndUnit("/accelerator/leavesB", this);
	leavesB->SetDefaultUnit("mm");
	leavesB->SetDefaultValue(300.);
	pAcc1->setLeavesBx(300.*mm);

	aperture1X=new G4UIcmdWithADoubleAndUnit("/accelerator/aperture1X", this);
	aperture1X->SetDefaultUnit("mm");
	aperture1X->SetDefaultValue(100.);
	pAcc1->setJaw1X(100.*mm);

	aperture1Y=new G4UIcmdWithADoubleAndUnit("/accelerator/aperture1Y", this);
	aperture1Y->SetDefaultUnit("mm");
	aperture1Y->SetDefaultValue(100.);
	pAcc1->setJaw1Y(100.*mm);

	aperture2X=new G4UIcmdWithADoubleAndUnit("/accelerator/aperture2X", this);
	aperture2X->SetDefaultUnit("mm");
	aperture2X->SetDefaultValue(-100.);
	pAcc1->setJaw2X(-100.*mm);

	aperture2Y=new G4UIcmdWithADoubleAndUnit("/accelerator/aperture2Y", this);
	aperture2Y->SetDefaultUnit("mm");
	aperture2Y->SetDefaultValue(-100.);
	pAcc1->setJaw2Y(-100.*mm);
}

CML2Acc1Messenger::~CML2Acc1Messenger(void)
{
	delete idEnergy;
	delete aperture1X;
	delete aperture2X;
	delete aperture1Y;
	delete aperture2Y;
	delete leavesA;
	delete leavesB;
}
void CML2Acc1Messenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
	if (cmd==aperture1X)
	{
		aperture1X->GetNewUnitValue(newValue);
		pAcc1->setJaw1X(aperture1X->GetNewDoubleValue(newValue));
	}

	if (cmd==aperture1Y)
	{
		aperture1Y->GetNewUnitValue(newValue);
		pAcc1->setJaw1Y(aperture1Y->GetNewDoubleValue(newValue));
	}

	if (cmd==aperture2X)
	{
		aperture2X->GetNewUnitValue(newValue);
		pAcc1->setJaw2X(aperture2X->GetNewDoubleValue(newValue));
	}
	if (cmd==aperture2Y)
	{
		aperture2Y->GetNewUnitValue(newValue);
		pAcc1->setJaw2Y(aperture2Y->GetNewDoubleValue(newValue));
	}

	if (cmd==leavesA)
	{pAcc1->setLeavesAx(leavesA->GetNewDoubleValue(newValue));}

	if (cmd==leavesB)
	{pAcc1->setLeavesBx(leavesA->GetNewDoubleValue(newValue));}


	if (cmd==idEnergy)
	{pAcc1->setidEnergy(idEnergy->GetNewIntValue(newValue));}
}
