#ifndef VDetectorOrganization_h
#define VDetectorOrganization_h

#include "G4Step.hh"

class G4StepPoint;

class VDetectorOrganization{

public:
  VDetectorOrganization(){};
  virtual ~VDetectorOrganization(){};   
  virtual unsigned int GetUnitID(const  G4Step* aStep) const =0;
  virtual int  Levels(const G4Step*) const;
  virtual void DetectorLevel(const G4Step*, int&, int*, G4String*) const;
};

#endif
