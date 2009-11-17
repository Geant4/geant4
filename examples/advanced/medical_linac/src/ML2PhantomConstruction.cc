#include "ML2PhantomConstruction.h"
#include "ML2PhantomConstructionMessenger.h"
#include "G4UImanager.hh"

CML2PhantomConstruction::CML2PhantomConstruction(void)
{
	this->phantomContstructionMessenger=new CML2PhantomConstructionMessenger(this);
}

CML2PhantomConstruction::~CML2PhantomConstruction(void)
{
	if (this->phantomName=="fullWater")
	{
		delete Ph_fullWater;
	}
	else if (this->phantomName=="BoxInBox")
	{
		delete Ph_BoxInBox;
	}
}

CML2PhantomConstruction* CML2PhantomConstruction::instance = 0;

CML2PhantomConstruction* CML2PhantomConstruction::GetInstance(void)
{
  if (instance == 0)
    {
      instance = new CML2PhantomConstruction();
     
    }
  return instance;
}
void CML2PhantomConstruction::design(void)
{
// switch between two different phantoms according to the main macro selection
	if (this->phantomName=="fullWater")
	{
		this->Ph_fullWater=new CML2Ph_FullWater();
	}
	else if (this->phantomName=="BoxInBox")
	{
		this->Ph_BoxInBox=new CML2Ph_BoxInBox();
	}

// read the messenger data related to the phantom selected 
	G4UImanager* UI = G4UImanager::GetUIpointer();
	G4String command = "/control/execute ";
	UI->ApplyCommand(command+this->PhantomSpecficationsFileName); // DEVE LEGGERE TUTTI I MESSENGER DELLE MACRO STRUTTURE
}
void CML2PhantomConstruction::Construct(G4VPhysicalVolume *PVWorld, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG)
{
// a call to select the right phantom
	this->design();
// create the phantom-world box
	G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4ThreeVector halfSize, centre;
	centre.set(0.*mm, 0.*mm, 0.*mm);
	halfSize.set(300.*mm, 300*mm, 300*mm);
	G4Box *phmWorldB = new G4Box("phmWorldG", halfSize.getX(), halfSize.getY(), halfSize.getZ());
	G4LogicalVolume *phmWorldLV = new G4LogicalVolume(phmWorldB, Vacuum, "phmWorldL", 0, 0, 0);
	G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::White());
	simpleAlSVisAtt->SetVisibility(true);
	simpleAlSVisAtt->SetForceWireframe(true);
	phmWorldLV->SetVisAttributes(simpleAlSVisAtt);

	G4RotationMatrix *rotation=new G4RotationMatrix();
	this->PVPhmWorld= new G4PVPlacement(rotation, centre, "phmWorldPV", phmWorldLV, PVWorld, false, 0);

// create the actual phantom
	if (this->phantomName=="fullWater")
	{
		this->Ph_fullWater->Construct(this->PVPhmWorld, saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG);
	}
	else if (this->phantomName=="BoxInBox")
	{
		this->Ph_BoxInBox->Construct(this->PVPhmWorld, saving_in_ROG_Voxels_every_events, seed, ROGOutFile, bSaveROG);
	}
}

