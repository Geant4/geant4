///////////////////////////////////////////////////////////////////////////////
// File: CCalGeometryConfiguration.cc
// Description: Handles geometry configuration to be constructed
///////////////////////////////////////////////////////////////////////////////

#include "CCalGeometryConfiguration.hh"

#include <fstream>

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
    cerr << "ERROR: In CCalGeometryConfiguration::getConstructFlag(const G4String& n)" 
	 << endl 
	 << "       " << n << " not found in configuration file" << endl;
  }

  return flag;
}

G4String CCalGeometryConfiguration::getFileName(const G4String& n) /*const*/ {
  G4String fn;
  CCalGeometryConfIterator it = theConfiguration.find(n);

  if (it != theConfiguration.end())
    fn = (*it).second.FileName;
  else {
    cerr << "ERROR: In CCalGeometryConfiguration::getConstructFlag(const G4String& n)" 
	 << endl 
	 << "       " << n << " not found in configuration file" << endl;
  }

  return fn;
}

CCalGeometryConfiguration::CCalGeometryConfiguration():
  theConfiguration() {

  ///////////////////////////////////////////////////////////////
  // Open the file
  G4String pathName = getenv("OSCAR_CONFPATH");
  G4String fileenv  = getenv("OSCAR_GEOMETRYCONF");
  if (!pathName || !fileenv) {
    cerr << "ERROR: OSCAR_GEOMETRYCONF and/or OSCAR_CONFPATH not set" << endl
	 << "       Set them to the geometry configuration file/path" << endl;
    exit(-2);
  }

  cout << " ==> Opening file " << fileenv << "..." << endl;
  ifstream is;
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
    cout << "CCalGeometryConfiguration constructor: Read \"" << name 
	 << "\" \"" << gcinfo.FileName << "\"" << tab << gcinfo.ConstructFlag 
	 << endl;
#endif
    theConfiguration[name] = gcinfo;
  }

  

  ///////////////////////////////////////////////////////////////
  // Close the file  
  is.close();
  cout << " <== Closed file " << fileenv << endl;
}
