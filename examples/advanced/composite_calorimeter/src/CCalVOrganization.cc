//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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
