///////////////////////////////////////////////////////////////////////////////
// File: CCalEcalOrganization.hh
// Description: Defines numbering schema for the Electromagnetic Calorimeter
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalEcalOrganization_h
#define CCalEcalOrganization_h

#include "CCalVOrganization.hh"

class CCalEcalOrganization: public CCalVOrganization {

public:
  CCalEcalOrganization(){};
  ~CCalEcalOrganization();
         
  virtual unsigned int GetUnitID(const G4Step* aStep) const ;
      
};

#endif
