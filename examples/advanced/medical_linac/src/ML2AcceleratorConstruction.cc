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


#include "ML2AcceleratorConstruction.hh"
#include "ML2AcceleratorConstructionMessenger.hh"


CML2AcceleratorConstruction::CML2AcceleratorConstruction(void)
{
	this->acceleratorConstructionMessenger=new CML2AcceleratorConstructionMessenger(this);
	this->idCurrentRotationX=0;
}

CML2AcceleratorConstruction::~CML2AcceleratorConstruction(void)
{
	if (this->AcceleratorName=="acc1")
	{delete this->accelerator1;}

	delete this->PVAccWorld;
	delete this->acceleratorConstructionMessenger;
}
CML2AcceleratorConstruction* CML2AcceleratorConstruction::instance = 0;

CML2AcceleratorConstruction* CML2AcceleratorConstruction::GetInstance(void)
{
  if (instance == 0)
    {
      instance = new CML2AcceleratorConstruction();
    }
  return instance;
}

void CML2AcceleratorConstruction::resetAccelerator()
{
	if (this->AcceleratorName=="acc1")
	{
		this->accelerator1->reset();
	}

}

bool CML2AcceleratorConstruction::design(void)
{
// switch between different accelerators according to the main macro selection (actually only one accelerator is available)
	std::cout << "I'm building "<< this->AcceleratorName<<"  accelerator"<< G4endl;
	bool bAccExists=false;
	if (this->AcceleratorName=="acc1")
	{this->accelerator1=CML2Acc1::GetInstance();bAccExists=true;}

	if (bAccExists && this->AcceleratorMacFileName!="")
	{
	// read the messenger data related to the accelerator selected 
		G4UImanager* UI = G4UImanager::GetUIpointer();
		G4String command = "/control/execute ";
		UI->ApplyCommand(command+this->AcceleratorMacFileName); 
	}

 	if (this->rotationsX.size()<1)
 	{this->addAcceleratorRotationsX(0.);}

	return bAccExists;
}
bool CML2AcceleratorConstruction::Construct(G4VPhysicalVolume *PVWorld, G4bool bOnlyVisio)
{
// a call to select the right accelerator
	this->bOnlyVisio=bOnlyVisio;
	if (this->design())
	{
		this->acceleratorConstructionMessenger->SetReferenceWorld(bOnlyVisio);
		// create the accelerator-world box
		G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
		G4ThreeVector halfSize;
		this->initialCentre.set(0.*mm, 0.*mm, -this->isoCentre);
		halfSize.set(600.*mm, 600.*mm, 600.*mm);
		G4Box *accWorldB = new G4Box("accWorldG", halfSize.getX(), halfSize.getY(), halfSize.getZ());
		G4LogicalVolume *accWorldLV = new G4LogicalVolume(accWorldB, Vacuum, "accWorldL", 0, 0, 0);
		G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::White());
		simpleAlSVisAtt->SetVisibility(false);
// 		simpleAlSVisAtt->SetForceWireframe(false);
		accWorldLV->SetVisAttributes(simpleAlSVisAtt);
	
		this->PVAccWorld= new G4PVPlacement(0, this->initialCentre, "acceleratorBox", accWorldLV, PVWorld, false, 0);

	// create the actual accelerator
		if (this->AcceleratorName=="acc1")
		{
			this->accelerator1->Construct(PVAccWorld, this->isoCentre);
			this->Z_Value_PhaseSpaceBeforeJaws=this->accelerator1->getBeforeJaws_Z_PhaseSpacePosition();
			this->accelerator1->writeInfo();
		}
	}
	else
	{
		return false;
	}
	return true;
}

void CML2AcceleratorConstruction::writeInfo()
{
	if (!this->bOnlyVisio)
	{std::cout <<"Actual rotation: "<<this->idCurrentRotationX<<"/"<<this->rotationsX.size() <<"  "<< G4endl;}
	std::cout <<"Accelerator angle: "<< this->currentRotationX/deg << " [deg]"<< G4endl;
}

G4RotationMatrix * CML2AcceleratorConstruction::rotateAccelerator()
{
	G4RotationMatrix *rmInv=new G4RotationMatrix();
	if (this->idCurrentRotationX <(int) this->rotationsX.size())
	{
		this->currentRotationX=this->rotationsX[this->idCurrentRotationX];
		rmInv=this->rotateAccelerator(this->currentRotationX);
		this->idCurrentRotationX++;
	}
	else
	{rmInv=0;}
	return rmInv;
}
G4RotationMatrix * CML2AcceleratorConstruction::rotateAccelerator(G4double angleX)
{
	this->currentRotationX=angleX;
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4ThreeVector NewCentre;
	G4RotationMatrix *rm=new G4RotationMatrix();
	G4RotationMatrix *rmInv=new G4RotationMatrix();
	this->PVAccWorld->SetTranslation(this->initialCentre);
	this->PVAccWorld->SetRotation(rm);
	if (this->bRotate90Y)	{rm->rotateY(90.*deg);}
	rm->rotateX(-angleX);
	this->PVAccWorld->SetRotation(rm);
	*rmInv=CLHEP::inverseOf(*rm);
	NewCentre=*rmInv*this->initialCentre;
	this->PVAccWorld->SetTranslation(NewCentre);
	G4GeometryManager::GetInstance()->CloseGeometry();
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
	return rmInv;
}
G4String CML2AcceleratorConstruction::getCurrentRotationString()
{
	char cR[5];
	G4int cRI=(G4int)((this->currentRotationX/deg)+.5);
	sprintf(cR,"%d",cRI);
	G4String rotationName=G4String(cR);
	if (this->bRotate90Y)
	{rotationName="_Ro90Y"+rotationName;}
	else
	{rotationName="_Ro"+rotationName;}
	return rotationName;
}
