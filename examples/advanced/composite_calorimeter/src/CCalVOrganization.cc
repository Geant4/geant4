///////////////////////////////////////////////////////////////////////////////
// File: CCalVOrganization.cc
// Description: Base class for definition of sensitive unit numbering schema
///////////////////////////////////////////////////////////////////////////////
#include "CCalVOrganization.hh"
#include "globals.hh"

int CCalVOrganization::Levels(const G4Step* aStep) const {

  //Find number of levels
  G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
  int level = 0;
  while (pv > 0) {
    level++;
    pv = pv->GetMother();
  }
  return level;
}


void CCalVOrganization::DetectorLevel(const G4Step* aStep, int& level,
				     int* copyno, G4String* name) const {

  //Get name and copy numbers
  if (level > 0) {
    G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
    for (int ii = level-1; ii >= 0; ii--) {
      name[ii]   = pv->GetName();
      copyno[ii] = pv->GetCopyNo();
      pv         = pv->GetMother();
    }
  }
}
