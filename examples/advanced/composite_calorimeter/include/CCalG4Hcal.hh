///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Hcal.hh
// Description: Euipped to construct the G4 geometry of the hadron calorimeter
//              in the 96 test beam run
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalG4Hcal_h
#define CCalG4Hcal_h 1

#include "CCalHcal.hh"
#include "CCalG4Able.hh"
#include "g4std/vector"

typedef G4LogicalVolume* ptrG4Log;

class CCalG4Hcal: public CCalHcal, public CCalG4Able {
public:
  //Constructor and Destructor
  CCalG4Hcal(const G4String &name);
  virtual ~CCalG4Hcal();

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
  G4std::vector<ptrG4Log> allSensitiveLogs;
};

#endif
