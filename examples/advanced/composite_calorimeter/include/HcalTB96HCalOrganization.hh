///////////////////////////////////////////////////////////////////////////////
// File: HcalTB96HCalOrganization.hh
// Date: 08/00
// Modifications:
///////////////////////////////////////////////////////////////////////////////
#ifndef HcalTB96HCalOrganization_h
#define HcalTB96HCalOrganization_h

#include "CaloOrganization.hh"

class HcalTB96HCalOrganization: public CaloOrganization{

public:
  HcalTB96HCalOrganization(){};
  ~HcalTB96HCalOrganization();
         
  virtual unsigned int GetUnitID(const G4Step* aStep) const ;
      
};

#endif
