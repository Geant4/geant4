///////////////////////////////////////////////////////////////////////////////
// File: CaloOrganization.hh
// Date: 29.10.99 
// Description: Definition of sensitive unit numbering schema
// Modifications:
///////////////////////////////////////////////////////////////////////////////
#ifndef CaloOrganization_h
#define CaloOrganization_h

#include "G4Step.hh"
#include "CMSCaloOrganization.hh"

class CaloOrganization {

public:
  CaloOrganization(){};
  virtual ~CaloOrganization(){};
	 
  virtual unsigned int GetUnitID(const G4Step* aStep) const = 0;
  virtual int  Levels(const G4Step*) const;
  virtual void DetectorLevel(const G4Step*, int&, int*, G4String*) const;
     
protected:
  CMSCaloOrganization theOrg;
};

#endif
