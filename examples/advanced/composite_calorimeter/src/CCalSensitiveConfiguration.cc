///////////////////////////////////////////////////////////////////////////////
// File: CCalSensitiveConfiguration.cc
// Description: CCalSensitiveConfiguration handles the declaration of
//              sensitive detectors and visibilities
///////////////////////////////////////////////////////////////////////////////
#include "CCalSensitiveConfiguration.hh"

#include <fstream.h>
#include <stdlib.h>


//Comment/Uncomment next line to hide/show debug information
//#define debug


CCalSensitiveConfiguration * CCalSensitiveConfiguration::instance = 0;

CCalSensitiveConfiguration* CCalSensitiveConfiguration::getInstance(){
  if (!instance) 
    instance = new CCalSensitiveConfiguration;
  return instance;
}


int CCalSensitiveConfiguration::getSensitiveFlag(const G4String& n) /*const*/ {
  int flag = -1;
  CCalSensitiveConfIterator it = theConfiguration.find(n);

  if (it != theConfiguration.end())
    flag = (*it).second.SensitiveFlag;
  else {
    cerr << "ERROR: In CCalSensitiveConfiguration::getConstructFlag(const "
	 << "G4String& n)" << endl 
	 << "       " << n << " not found in configuration file" << endl;
  }

  return flag;
}

G4String CCalSensitiveConfiguration::getFileName(const G4String& n) /*const*/ {
  G4String fn;
  CCalSensitiveConfIterator it = theConfiguration.find(n);

  if (it != theConfiguration.end())
    fn = (*it).second.FileName;
  else {
    cerr << "ERROR: In CCalSensitiveConfiguration::getConstructFlag(const "
	 << "G4String& n)" << endl 
	 << "       " << n << " not found in configuration file" << endl;
  }

  return fn;
}

CCalSensitiveConfiguration::CCalSensitiveConfiguration():
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
    cout << "CCalSensitiveConfiguration constructor: Read \"" << name 
	 << "\" \"" << gcinfo.FileName << "\"" << tab << gcinfo.SensitiveFlag 
	 << endl;
#endif
    theConfiguration[name] = gcinfo;
  }
  

  ///////////////////////////////////////////////////////////////
  // Close the file  
  is.close();
  cout << " <== Closed file " << fileenv << endl;
}
