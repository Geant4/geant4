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
#include "CCalEcalOrganization.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "SystemOfUnits.h"

CCalEcalOrganization::~CCalEcalOrganization(){
  G4cout << " Deleting CCalEcalOrganization" << G4endl;
}


unsigned int CCalEcalOrganization::GetUnitID(const G4Step* aStep) const {

  G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
  int idx=0, idl=0, idn=0;
  if (pv > 0) {
    idx = pv->GetCopyNo() - 1;
    pv  = pv->GetMother();
    if (pv > 0) {
      idl = pv->GetCopyNo() - 1;
      pv  = pv->GetMother();
      if (pv > 0) 
	idn = pv->GetCopyNo();
    }
  }
  unsigned int idunit = idn*4096 + idl*64 + idx;
  
  return idunit;
}
