///////////////////////////////////////////////////////////////////////////////
// File: CCalVOrganization.hh
// Description: Base class for definition of sensitive unit numbering schema
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalVOrganization_h
#define CCalVOrganization_h

#include "G4Step.hh"
#include "CCaloOrganization.hh"

class CCalVOrganization {

public:
  CCalVOrganization(){};
  virtual ~CCalVOrganization(){};
	 
  virtual unsigned int GetUnitID(const G4Step* aStep) const = 0;
  virtual int  Levels(const G4Step*) const;
  virtual void DetectorLevel(const G4Step*, int&, int*, G4String*) const;
     
protected:
  CCaloOrganization theOrg;
};

#endif
