///////////////////////////////////////////////////////////////////////////////
// File: CrystalMatrixOrganization.hh
// Date: 08/00 
// Modifications:
///////////////////////////////////////////////////////////////////////////////
#ifndef CrystalMatrixOrganization_h
#define CrystalMatrixOrganization_h

#include "CaloOrganization.hh"

class CrystalMatrixOrganization: public CaloOrganization {

public:
  CrystalMatrixOrganization(){};
  ~CrystalMatrixOrganization();
         
  virtual unsigned int GetUnitID(const G4Step* aStep) const ;
      
};

#endif
