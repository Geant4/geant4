#include "ML2AcceleratorConstructionMessenger.h"
#include "ML2AcceleratorConstruction.h"
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
