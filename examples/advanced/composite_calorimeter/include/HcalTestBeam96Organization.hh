#ifndef HcalTestBeam96Organization_h
#define HcalTestBeam96Organization_h

#include "VDetectorOrganization.hh"

class G4Step;
class HcalTestBeam96Organization: public VDetectorOrganization{

public:
  HcalTestBeam96Organization();
  ~HcalTestBeam96Organization();
	 
  virtual unsigned int GetUnitID(const G4Step* aStep) const ;

public:

};

#endif
