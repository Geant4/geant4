///////////////////////////////////////////////////////////////////////////////
// File: G4GeometryConfiguration.h
// Date: 03/98 
// Modifications: 24/03/98 I.G.
//                03/05/99 I.G. Added sensitive behaviour
//                16/12/99 I.G. Updated for STL.  Needs revision!!! 
//                   03/00 S.B. in OSCAR
// Description: This singleton holds the information given in the file 
//              g4geometry.conf 
//              Use getInstance to retrieve it.
///////////////////////////////////////////////////////////////////////////////
#ifndef G4GeometryConfiguration_h
#define G4GeometryConfiguration_h 1

//#include <rw/tvhdict.h>
#include "g4std/map"
#include "utils.hh"

class G4GeometryConfiguration {

  struct GCInfo { 
    G4String FileName;
    int      SensitiveFlag;
  };

   typedef G4std::map<G4String, GCInfo, less<G4String> > G4GeometryConfTable;
   typedef G4std::map<G4String, GCInfo, less<G4String> >::iterator G4GeometryConfIterator;
		    
  
public:
  ~G4GeometryConfiguration() {}

  static G4GeometryConfiguration* getInstance();

  int      getSensitiveFlag(const G4String& n) /*const*/;
  G4String getFileName(const G4String& n) /*const*/;

private:
  G4GeometryConfiguration();

private:
  static G4GeometryConfiguration* instance;
  G4GeometryConfTable theConfiguration;
};

#endif
