///////////////////////////////////////////////////////////////////////////////
// File: CCalHcalOrganization.hh
// Description: Defines numbering schema for the Hadron Calorimeter
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalHcalOrganization_h
#define CCalHcalOrganization_h

#include "CCalVOrganization.hh"

class CCalHcalOrganization: public CCalVOrganization{

public:
  CCalHcalOrganization(){};
  ~CCalHcalOrganization();
         
  virtual unsigned int GetUnitID(const G4Step* aStep) const ;
      
};

#endif
