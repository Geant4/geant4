///////////////////////////////////////////////////////////////////////////////
// File: CaloOrganization.hh
// Date: 29.10.99 
// Description: Definition of sensitive unit numbering schema
// Modifications:
///////////////////////////////////////////////////////////////////////////////
#ifndef CaloOrganization_h
#define CaloOrganization_h

#include "VDetectorOrganization.hh"
#include "CMSCaloOrganization.hh"

class CaloOrganization: public VDetectorOrganization{

public:
  CaloOrganization(){};
  virtual ~CaloOrganization(){};
	 
  virtual unsigned int GetUnitID(const G4Step* aStep) const {return 0;}
      
protected:
  CMSCaloOrganization theOrg;
};

#endif
