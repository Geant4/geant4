///////////////////////////////////////////////////////////////////////////////
// File: G4HcalTB96HCal.hh
// Date: 08/0
// Modifications: 
// Description: Euipped to construct the G4 geometry of the hadron calorimeter
//              in the 96 test beam run
///////////////////////////////////////////////////////////////////////////////
#ifndef G4HcalTB96HCal_h
#define G4HcalTB96HCal_h 1

#include "HcalTB96HCal.hh"
#include "CCalG4Able.hh"
#include <vector>

typedef G4LogicalVolume* ptrG4Log;

class G4HcalTB96HCal: public HcalTB96HCal, public CCalG4Able {
public:
  //Constructor and Destructor
  G4HcalTB96HCal(const G4String &name);
  virtual ~G4HcalTB96HCal();

protected:
  //This methods actually constructs the volume.
  virtual G4VPhysicalVolume* constructIn(G4VPhysicalVolume*);
  virtual void constructDaughters();

  //Construct layer of scintillator or absorber
  G4LogicalVolume*   constructScintillatorLayer (G4int);
  G4LogicalVolume*   constructAbsorberLayer (G4int);

  //Constructs the sensitive detectors and associates them to the corresponding
  //logical volumes
  virtual void constructSensitive() ;

private:
  //Private data members
  ptrG4Log*        sclLog;
  ptrG4Log*        absLog;

  // Logical volumes for sensitive detectors
  vector<ptrG4Log> allSensitiveLogs;
};

#endif
