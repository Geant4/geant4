///////////////////////////////////////////////////////////////////////////////
// File: CCalHcalOrganization.cc
// Description: Defines numbering schema for the Hadron Calorimeter
///////////////////////////////////////////////////////////////////////////////

#include "CCalHcalOrganization.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "SystemOfUnits.h"

CCalHcalOrganization::~CCalHcalOrganization() {
  G4cout << " Deleting CCalHcalOrganization" << G4endl;
}


unsigned int CCalHcalOrganization::GetUnitID(const G4Step* aStep) const {

  G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
  if (pv > 0) 
    pv = pv->GetMother();
  int idunit=0;
  if (pv > 0)
    idunit = pv->GetCopyNo();

  return idunit;
}
