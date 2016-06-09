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


#include "ML2AcceleratorConstruction.hh"
#include "ML2AcceleratorConstructionMessenger.hh"
#include "G4UImanager.hh"


CML2AcceleratorConstruction::CML2AcceleratorConstruction(void)
{
	this->acceleratorConstructionMessenger=new CML2AcceleratorConstructionMessenger(this);
}

CML2AcceleratorConstruction::~CML2AcceleratorConstruction(void)
{
	delete this->PVAccWorld;
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
void CML2AcceleratorConstruction::design(void)
{
// switch between two different phantoms according to the main macro selection (actually only one accelerator is available

	if (this->AcceleratorName=="acc1")
	{this->accelerator1=CML2Acc1::GetInstance();}

// read the messenger data related to the accelerator selected 
	G4UImanager* UI = G4UImanager::GetUIpointer();
	G4String command = "/control/execute ";
	UI->ApplyCommand(command+this->AcceleratorSpecficationsFileName); 
}

void CML2AcceleratorConstruction::Construct(G4VPhysicalVolume *PVWorld)
{
// a call to select the right accelerator
	this->design();
// create the accelerator-world box
	G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4ThreeVector halfSize, centre;
	centre.set(0.*mm, 0.*mm, -1300.*mm);
	halfSize.set(800.*mm, 800.*mm, 1000.*mm);
	G4Box *accWorldB = new G4Box("accWorldG", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	G4LogicalVolume *accWorldLV = new G4LogicalVolume(accWorldB, Vacuum, "accWorldL", 0, 0, 0);
	G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::White());
	simpleAlSVisAtt->SetVisibility(true);
	simpleAlSVisAtt->SetForceWireframe(true);
	accWorldLV->SetVisAttributes(simpleAlSVisAtt);

	G4RotationMatrix *rotation=new G4RotationMatrix();
	centre=*rotation*centre;
	centre.setY(-centre.getY());

	this->PVAccWorld= new G4PVPlacement(rotation, centre, "acceleratorBox", accWorldLV, PVWorld, false, 0);

// create the actual accelerator
	if (this->AcceleratorName=="acc1")
	{this->accelerator1->Construct(PVAccWorld);}
}

