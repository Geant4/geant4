///////////////////////////////////////////////////////////////////////////////
// File: CMSSensitiveDetectors.h
// Date: 15/05/02 S.B.
// Modifications: 
// Description: Container of logical volume pointers which can be sensitive
//              detectors
///////////////////////////////////////////////////////////////////////////////
#ifndef CMSSensitiveDetectors_h
#define CMSSensitiveDetectors_h 1

#include <vector>
#include <map>
#include "G4VSensitiveDetector.hh"
#include "G4LogicalVolume.hh"

typedef multimap< string, G4LogicalVolume*, less<string> > mmslv;

class CMSSensitiveDetectors {

public:    
  ~CMSSensitiveDetectors(){};
  vector<G4LogicalVolume*> getVolumes (const string& name, bool exists = 0);
  void registerVolume (const string& name, G4LogicalVolume*);
  bool setSensitive(const string& string, G4VSensitiveDetector* sens);
  static CMSSensitiveDetectors* getInstance();
private:
  CMSSensitiveDetectors(){};
  static CMSSensitiveDetectors* theInstance;
  // logical volume container: name, G4LogicalVolume*
  mmslv theLVs;

};

#endif
