///////////////////////////////////////////////////////////////////////////////
// File: CCalSensAssign.hh
// Description: CCalSenAssign creates and assigns the sensitive detetctors 
//              from the map of logical volumes which are potentially sensitive
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSensAssign_h
#define CCalSensAssign_h

#include "CaloOrganization.hh"

#include <map>
#include "globals.hh"
#include "G4VSensitiveDetector.hh"

class CCalSensAssign {

public:
  ~CCalSensAssign(){};
  static CCalSensAssign* getInstance();
  bool assign();
  bool stackingAction();
  bool addCaloSD(G4String name, CaloOrganization* numberingScheme);

private:
  CCalSensAssign();

  static CCalSensAssign* theInstance;
  map<G4String,G4VSensitiveDetector*> sens_;

};
#endif
