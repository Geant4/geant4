///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Hall.hh
// Description: Equipped to construct G4 geometry of the experimental hall
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalG4Hall_h
#define CCalG4Hall_h 1

#include "CCalHall.hh"
#include "CCalG4Able.hh"

class CCalG4Hall: public CCalHall, public CCalG4Able {
public:
  //Constructor and Destructor
  CCalG4Hall(const G4String &name);
  virtual ~CCalG4Hall();
  
protected:
  //This methods actually constructs the volume.
  virtual G4VPhysicalVolume* constructIn(G4VPhysicalVolume*);
  virtual void constructDaughters();
};

#endif
