///////////////////////////////////////////////////////////////////////////////
// File: G4GeometryConfiguration.cc
// Date: 03/98 I. Gonzalez
// Modifications: 27/03/00 SB In OSCAR
//                02/05/00 I.G. S.B. Use local file if available
///////////////////////////////////////////////////////////////////////////////
#include "G4GeometryConfiguration.hh"

#include <fstream.h>
#include <stdlib.h>


//Comment/Uncomment next line to hide/show debug information
//#define debug


G4GeometryConfiguration * G4GeometryConfiguration::instance = 0;

G4GeometryConfiguration* G4GeometryConfiguration::getInstance(){
  if (!instance) 
    instance = new G4GeometryConfiguration;
  return instance;
}


int G4GeometryConfiguration::getSensitiveFlag(const G4String& n) /*const*/ {
  int flag = -1;
  G4GeometryConfIterator it = theConfiguration.find(n);

  if (it != theConfiguration.end())
    flag = (*it).second.SensitiveFlag;
  else {
    cerr << "ERROR: In GeometryConfiguration::getConstructFlag(const G4String& n)" 
	 << endl 
	 << "       " << n << " not found in configuration file" << endl;
  }

  return flag;
}

G4String G4GeometryConfiguration::getFileName(const G4String& n) /*const*/ {
  G4String fn;
  G4GeometryConfIterator it = theConfiguration.find(n);

  if (it != theConfiguration.end())
    fn = (*it).second.FileName;
  else {
    cerr << "ERROR: In GeometryConfiguration::getConstructFlag(const G4String& n)" 
	 << endl 
	 << "       " << n << " not found in configuration file" << endl;
  }

  return fn;
}

G4GeometryConfiguration::G4GeometryConfiguration():
  theConfiguration() {

  ///////////////////////////////////////////////////////////////
  // Open the file
  G4String pathName = getenv("OSCAR_CONFPATH");
  G4String fileenv  = getenv("OSCAR_SENSITIVECONF");
  if (!pathName || !fileenv) {
    cerr << "ERROR: OSCAR_SENSITIVECONF and/or OSCAR_CONFPATH not set" << endl
	 << "       Set them to the sensitive configuration file/path" << endl;
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
    is >> gcinfo.SensitiveFlag >> jump;
#ifdef debug
    cout << "G4GeometryConfiguration constructor: Read \"" << name << "\" \"" 
	 << gcinfo.FileName << "\"" << tab << gcinfo.SensitiveFlag << endl;
#endif
    theConfiguration[name] = gcinfo;
  }

  

  ///////////////////////////////////////////////////////////////
  // Close the file  
  is.close();
  cout << " <== Closed file " << fileenv << endl;
}
