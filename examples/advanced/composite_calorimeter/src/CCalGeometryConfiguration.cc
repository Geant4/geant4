///////////////////////////////////////////////////////////////////////////////
// File: CCalGeometryConfiguration.cc
// Description: Handles geometry configuration to be constructed
///////////////////////////////////////////////////////////////////////////////

#include "CCalGeometryConfiguration.hh"

#include "g4std/fstream"

//Comment/Uncomment next line to hide/show debug information
//#define debug


CCalGeometryConfiguration * CCalGeometryConfiguration::instance = 0;

CCalGeometryConfiguration* CCalGeometryConfiguration::getInstance(){
  if (!instance) 
    instance = new CCalGeometryConfiguration;
  return instance;
}


int CCalGeometryConfiguration::getConstructFlag(const G4String& n) /*const*/ {
  int flag = -1;
  CCalGeometryConfIterator it = theConfiguration.find(n);

  if (it != theConfiguration.end())
    flag = (*it).second.ConstructFlag;
  else {
    G4cerr << "ERROR: In CCalGeometryConfiguration::getConstructFlag(const G4String& n)" 
	 << G4endl 
	 << "       " << n << " not found in configuration file" << G4endl;
  }

  return flag;
}

G4String CCalGeometryConfiguration::getFileName(const G4String& n) /*const*/ {
  G4String fn;
  CCalGeometryConfIterator it = theConfiguration.find(n);

  if (it != theConfiguration.end())
    fn = (*it).second.FileName;
  else {
    G4cerr << "ERROR: In CCalGeometryConfiguration::getConstructFlag(const G4String& n)" 
	 << G4endl 
	 << "       " << n << " not found in configuration file" << G4endl;
  }

  return fn;
}

CCalGeometryConfiguration::CCalGeometryConfiguration():
  theConfiguration() {

  ///////////////////////////////////////////////////////////////
  // Open the file
  G4String pathName = getenv("CCAL_CONFPATH");
  G4String fileenv  = getenv("CCAL_GEOMETRYCONF");
  if (!pathName || !fileenv) {
    G4cerr << "ERROR: CCAL_GEOMETRYCONF and/or CCAL_CONFPATH not set" << G4endl
	 << "       Set them to the geometry configuration file/path" << G4endl;
    exit(-2);
  }

  G4cout << " ==> Opening file " << fileenv << "..." << G4endl;
  G4std::ifstream is;
  bool ok = openGeomFile(is, pathName, fileenv);
  if (!ok)
    exit(-1);

  G4String name;
  GCInfo gcinfo;

  while (is) {
    readName(is, name);
    readName(is, gcinfo.FileName);
    is >> gcinfo.ConstructFlag >> jump;
#ifdef debug
    G4cout << "CCalGeometryConfiguration constructor: Read \"" << name 
	 << "\" \"" << gcinfo.FileName << "\"" << tab << gcinfo.ConstructFlag 
	 << G4endl;
#endif
    theConfiguration[name] = gcinfo;
  }

  

  ///////////////////////////////////////////////////////////////
  // Close the file  
  is.close();
  G4cout << " <== Closed file " << fileenv << G4endl;
}
