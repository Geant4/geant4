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
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

CML2AcceleratorConstruction::CML2AcceleratorConstruction(void)
{
	acceleratorConstructionMessenger = new CML2AcceleratorConstructionMessenger(this);
	idCurrentRotationX = 0;
}

CML2AcceleratorConstruction::~CML2AcceleratorConstruction(void)
{
	delete acceleratorConstructionMessenger;
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
		accelerator -> reset();
}

bool CML2AcceleratorConstruction::design(void)
{
// switch between different accelerators according to the main macro selection (actually only one accelerator is available)
	G4cout << "I'm building " << AcceleratorName << " accelerator" << G4endl;
	bool bAccExists = false;
	if (AcceleratorName == "acc1")
	{
		accelerator = CML2Acc1::GetInstance();
		bAccExists = true;
	}
	else if (AcceleratorName == "acc2")
	{
		accelerator = CML2Acc2::GetInstance();
		bAccExists = true;
	}
	else if (AcceleratorName == "accSaturn")
	{
		accelerator = CML2AccSaturn::GetInstance();
		bAccExists = true;
	}

	if (bAccExists && AcceleratorMacFileName!="")
	{
		// read the messenger data related to the selected accelerator
		G4UImanager* UI = G4UImanager::GetUIpointer();
		G4String command = "/control/execute ";
		UI->ApplyCommand(command+AcceleratorMacFileName); 
	}

 	if (rotationsX.size() < 1)
 	{
 		addAcceleratorRotationsX(0.);
 	}

	return bAccExists;
}
bool CML2AcceleratorConstruction::Construct(G4VPhysicalVolume *PVWorld, G4bool bOV)
{
	// a call to select the right accelerator
	bOnlyVisio = bOV;
	if (design())
	{
//		G4cout << "*** debug *** AcceleratorConstruction::Construct" << G4endl;
		acceleratorConstructionMessenger->SetReferenceWorld(bOnlyVisio);
		// create the accelerator-world box
		G4Material *Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
		G4ThreeVector halfSize;
		initialCentre.set(0.*mm, 0.*mm, -isoCentre);
		halfSize.set(600.*mm, 600.*mm, 600.*mm);
		G4Box *accWorldB = new G4Box("accWorldG", halfSize.getX(), halfSize.getY(), halfSize.getZ());
		G4LogicalVolume *accWorldLV = new G4LogicalVolume(accWorldB, Vacuum, "accWorldL", 0, 0, 0);
		G4VisAttributes* simpleAlSVisAtt = new G4VisAttributes(G4Colour::White());
		simpleAlSVisAtt -> SetVisibility(false);
		accWorldLV -> SetVisAttributes(simpleAlSVisAtt);
	
		PVAccWorld= new G4PVPlacement(0, initialCentre, "acceleratorBox", accWorldLV, PVWorld, false, 0);

		// create the actual accelerator
		accelerator -> Construct(PVAccWorld, isoCentre);
		Z_Value_PhaseSpaceBeforeJaws = accelerator -> getBeforeJaws_Z_PhaseSpacePosition();
		accelerator -> writeInfo();

	}
	else
	{
		return false;
	}
	return true;
}

void CML2AcceleratorConstruction::writeInfo()
{
	if (!bOnlyVisio)
	{
		G4cout << "Actual rotation: " << idCurrentRotationX << "/" << rotationsX.size() << "  " << G4endl;
	}
	G4cout << "Accelerator angle: " << currentRotationX/deg << " [deg]" << G4endl;
}

G4RotationMatrix * CML2AcceleratorConstruction::rotateAccelerator()
{
	G4RotationMatrix *rmInv=new G4RotationMatrix();
	if (idCurrentRotationX <(int) rotationsX.size())
	{
		currentRotationX = rotationsX[idCurrentRotationX];
		rmInv = rotateAccelerator(currentRotationX);
		idCurrentRotationX++;
	}
	else
	{
		rmInv = 0;
	}
	return rmInv;
}
G4RotationMatrix * CML2AcceleratorConstruction::rotateAccelerator(G4double angleX)
{
	currentRotationX = angleX;
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4ThreeVector NewCentre;
	G4RotationMatrix *rm = new G4RotationMatrix();
	G4RotationMatrix *rmInv = new G4RotationMatrix();
	PVAccWorld->SetTranslation(initialCentre);
	PVAccWorld->SetRotation(rm);
	if (bRotate90Y)
	{
		rm->rotateY(90.*deg);
	}
	rm->rotateX(-angleX);
	PVAccWorld->SetRotation(rm);
	*rmInv=CLHEP::inverseOf(*rm);
	NewCentre=*rmInv*initialCentre;
	PVAccWorld->SetTranslation(NewCentre);
	G4GeometryManager::GetInstance()->CloseGeometry();
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
	return rmInv;
}
G4String CML2AcceleratorConstruction::getCurrentRotationString()
{
	char cR[5];
	G4int cRI=(G4int)((currentRotationX/deg)+.5);
	sprintf(cR,"%d",cRI);
	G4String rotationName=G4String(cR);
	if (bRotate90Y)
	{rotationName="_Ro90Y"+rotationName;}
	else
	{rotationName="_Ro"+rotationName;}
	return rotationName;
}
