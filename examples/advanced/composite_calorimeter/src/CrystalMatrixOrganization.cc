///////////////////////////////////////////////////////////////////////////////
// File: CrystalMatrixOrganization.cc
// Date: 08/00 S.Banerjee
// Modifications:
///////////////////////////////////////////////////////////////////////////////
#include "CrystalMatrixOrganization.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "SystemOfUnits.h"

CrystalMatrixOrganization::~CrystalMatrixOrganization(){
  cout << " Deleting CrystalMatrixOrganization" << endl;
}


unsigned int CrystalMatrixOrganization::GetUnitID(const G4Step* aStep) const {

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
