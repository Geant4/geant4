///////////////////////////////////////////////////////////////////////////////
// File: HcalTB96HCalOrganization.cc
// Date: 08/00 S.Banerjee
// Modifications:
///////////////////////////////////////////////////////////////////////////////

#include "HcalTB96HCalOrganization.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "SystemOfUnits.h"

HcalTB96HCalOrganization::~HcalTB96HCalOrganization() {
  cout << " Deleting HcalTB96HCalOrganization" << endl;
}


unsigned int HcalTB96HCalOrganization::GetUnitID(const G4Step* aStep) const {

  G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
  if (pv > 0) 
    pv = pv->GetMother();
  int idunit=0;
  if (pv > 0)
    idunit = pv->GetCopyNo();

  return idunit;
}
