///////////////////////////////////////////////////////////////////////////////
// File: VDetectorOrganization.cc
// Author: Veronique Lefebure
// Modifications: 15/05/02 SB Add level of volume names and copy numbers
///////////////////////////////////////////////////////////////////////////////
#include "VDetectorOrganization.hh"
#include "globals.hh"

int VDetectorOrganization::Levels(const G4Step* aStep) const {

  //Find number of levels
  G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
  int level = 0;
  while (pv > 0) {
    level++;
    pv = pv->GetMother();
  }
  return level;
}


void VDetectorOrganization::DetectorLevel(const G4Step* aStep, int& level,
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
