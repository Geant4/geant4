///////////////////////////////////////////////////////////////////////////////
// File: CCalSensitiveDetectors.hh
// Description: Container of logical volume pointers which can be sensitive
//              detectors
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSensitiveDetectors_h
#define CCalSensitiveDetectors_h 1

#include <vector>
#include <map>
#include "G4VSensitiveDetector.hh"
#include "G4LogicalVolume.hh"

typedef multimap< string, G4LogicalVolume*, less<string> > mmslv;

class CCalSensitiveDetectors {

public:    
  ~CCalSensitiveDetectors(){};
  vector<G4LogicalVolume*> getVolumes (const string& name, bool exists = 0);
  void registerVolume (const string& name, G4LogicalVolume*);
  bool setSensitive(const string& string, G4VSensitiveDetector* sens);
  static CCalSensitiveDetectors* getInstance();
private:
  CCalSensitiveDetectors(){};
  static CCalSensitiveDetectors* theInstance;
  // logical volume container: name, G4LogicalVolume*
  mmslv theLVs;

};

#endif
