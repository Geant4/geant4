//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
