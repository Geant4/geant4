///////////////////////////////////////////////////////////////////////////////
// File: CCalSensitiveConfiguration.h
// Description: This singleton holds the information given in the file 
//              g4geometry.conf 
//              Use getInstance to retrieve it.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSensitiveConfiguration_h
#define CCalSensitiveConfiguration_h 1

#include "g4std/map"
#include "utils.hh"

class CCalSensitiveConfiguration {

  struct GCInfo { 
    G4String FileName;
    int      SensitiveFlag;
  };

   typedef G4std::map<G4String, GCInfo, less<G4String> > CCalSensitiveConfTable;
   typedef G4std::map<G4String, GCInfo, less<G4String> >::iterator CCalSensitiveConfIterator;
		    
  
public:
  ~CCalSensitiveConfiguration() {}

  static CCalSensitiveConfiguration* getInstance();

  int      getSensitiveFlag(const G4String& n) /*const*/;
  G4String getFileName(const G4String& n) /*const*/;

private:
  CCalSensitiveConfiguration();

private:
  static CCalSensitiveConfiguration* instance;
  CCalSensitiveConfTable theConfiguration;
};

#endif
