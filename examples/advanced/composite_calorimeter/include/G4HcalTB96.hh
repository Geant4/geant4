///////////////////////////////////////////////////////////////////////////////
// File: G4HcalTB96.hh
// Date: 08/00
// Modifications: 
///////////////////////////////////////////////////////////////////////////////
#ifndef G4HcalTB96_h
#define G4HcalTB96_h 1

#include "HcalTB96.hh"
#include "CCalG4Able.hh"

class G4HcalTB96: public HcalTB96, public CCalG4Able {
public:
  //Constructor and Destructor
  G4HcalTB96(const G4String &name);
  virtual ~G4HcalTB96();
  
protected:
  //This methods actually constructs the volume.
  virtual G4VPhysicalVolume* constructIn(G4VPhysicalVolume*);
  virtual void constructDaughters();
};

#endif
