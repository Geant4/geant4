#include "HcalTestBeam96Organization.hh"
#include "globals.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"

HcalTestBeam96Organization::HcalTestBeam96Organization(){}

HcalTestBeam96Organization::~HcalTestBeam96Organization(){

 cout <<"Deleting HcalTestBeam96Organization"<<endl;

}
	 
unsigned int HcalTestBeam96Organization::GetUnitID(const G4Step* aStep) const{

  G4StepPoint* PreStepPoint = aStep->GetPreStepPoint(); 
  int CopyNo   = PreStepPoint->GetPhysicalVolume()->GetCopyNo();
  return CopyNo;

}

