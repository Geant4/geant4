#ifndef CMSSensAssign_h
#define CMSSensAssign_h

#include "VDetectorOrganization.hh"

#include <map>
#include "globals.hh"
#include "G4VSensitiveDetector.hh"

class CMSSensAssign {

public:
  ~CMSSensAssign(){};
  static CMSSensAssign* getInstance();
  bool assign();
  bool stackingAction();
  bool addCaloSD(G4String name, VDetectorOrganization* numberingScheme);

private:
  CMSSensAssign();

  static CMSSensAssign* theInstance;
  map<G4String,G4VSensitiveDetector*> sens_;

};
#endif
