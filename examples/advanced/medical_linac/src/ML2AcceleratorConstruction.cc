#include "ML2AcceleratorConstruction.h"
#include "ML2AcceleratorConstructionMessenger.h"
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

