///////////////////////////////////////////////////////////////////////////////
// File: GeometryConfiguration.hh
// Date: 03/98 
// Modifications: 24/03/98 I.G.
//                03/05/99 I.G. Added sensitive behaviour
//                16/12/99 I.G. Updated for STL.  Needs revision!!! 
// Description: This singleton holds the information in the file geometry.conf 
//              Use getInstance to retrieve it.
///////////////////////////////////////////////////////////////////////////////
#ifndef GeometryConfiguration_h
#define GeometryConfiguration_h 1

#include "g4std/map"
#include "utils.hh"

class GeometryConfiguration {

  struct GCInfo { 
    G4String FileName;
    int      ConstructFlag;
  };

//  typedef RWTValHashDictionary<G4String,GCInfo> GeometryConfTable;
//  typedef RWTValHashDictionaryIterator<G4String,GCInfo> GeometryConfIterator;
  typedef G4std::map<G4String, GCInfo, less<G4String> > GeometryConfTable;
  typedef G4std::map<G4String, GCInfo, less<G4String> >::iterator GeometryConfIterator;
		    
  
public:
  ~GeometryConfiguration() {}

  static GeometryConfiguration* getInstance();

  int      getConstructFlag(const G4String& n) /*const*/;
  G4String getFileName(const G4String& n) /*const*/;

  GeometryConfTable& getConfTable() { 
    return theConfiguration;
  }

private:
  GeometryConfiguration();

private:
  static GeometryConfiguration* instance;
  GeometryConfTable theConfiguration;
};

#endif
