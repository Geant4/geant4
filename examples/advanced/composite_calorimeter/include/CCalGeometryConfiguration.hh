///////////////////////////////////////////////////////////////////////////////
// File: CCalGeometryConfiguration.hh
// Description: This singleton holds the information in the file geometry.conf 
//              Use getInstance to retrieve it.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalGeometryConfiguration_h
#define CCalGeometryConfiguration_h 1

#include "g4std/map"
#include "utils.hh"

class CCalGeometryConfiguration {

  struct GCInfo { 
    G4String FileName;
    int      ConstructFlag;
  };

  typedef G4std::map<G4String, GCInfo, less<G4String> > CCalGeometryConfTable;
  typedef G4std::map<G4String, GCInfo, less<G4String> >::iterator CCalGeometryConfIterator;
		    
  
public:
  ~CCalGeometryConfiguration() {}

  static CCalGeometryConfiguration* getInstance();

  int      getConstructFlag(const G4String& n) /*const*/;
  G4String getFileName(const G4String& n) /*const*/;

  CCalGeometryConfTable& getConfTable() { 
    return theConfiguration;
  }

private:
  CCalGeometryConfiguration();

private:
  static CCalGeometryConfiguration* instance;
  CCalGeometryConfTable theConfiguration;
};

#endif
