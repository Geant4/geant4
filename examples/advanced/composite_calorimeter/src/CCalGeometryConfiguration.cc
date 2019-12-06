//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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


G4int CCalGeometryConfiguration::getConstructFlag(const G4String& n) /*const*/ {
  G4int flag = -1;
  auto it = theConfiguration.find(n);

  if (it != theConfiguration.cend())
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
  auto it = theConfiguration.find(n);

  if (it != theConfiguration.cend())
    fn = (*it).second.FileName;
  else {
    G4cerr << "ERROR: In CCalGeometryConfiguration::getFileName(const G4String& n)" 
         << G4endl 
         << "       " << n << " not found in configuration file" << G4endl;
  }

  return fn;
}

CCalGeometryConfiguration::CCalGeometryConfiguration():
  theConfiguration() {

  ///////////////////////////////////////////////////////////////
  // Open the file
  G4String pathName = std::getenv("CCAL_CONFPATH");
  G4String fileenv  = std::getenv("CCAL_GEOMETRYCONF");
  if (!pathName || !fileenv) {
     G4ExceptionDescription ed;
     ed << "ERROR: CCAL_GEOMETRYCONF and/or CCAL_CONFPATH not set" << G4endl
        << "       Set them to the geometry configuration file/path" << G4endl;
     G4Exception("CCalGeometryConfiguration::CCalGeometryConfiguration()",
                 "ccal003",
                 FatalException,ed);
  }

  G4cout << " ==> Opening file " << fileenv << "..." << G4endl;
  std::ifstream is;
  G4bool ok = openGeomFile(is, pathName, fileenv);
  if (!ok)
    {
      G4Exception("CCalGeometryConfiguration::CCalGeometryConfiguration()",
                 "ccal004",
                 FatalException,"Unable to open input data file");
    }


  G4String name;
  GCInfo gcinfo;

  while (is) {
    readName(is, name);
    readName(is, gcinfo.FileName);
    is >> gcinfo.ConstructFlag >> jump;
#ifdef debug
    G4cout << "CCalGeometryConfiguration constructor: Read \"" << name 
         << " \"" << gcinfo.FileName << "\"" << tab << gcinfo.ConstructFlag 
         << G4endl;
#endif
    theConfiguration[name] = gcinfo;
  }

  

  ///////////////////////////////////////////////////////////////
  // Close the file  
  is.close();
  G4cout << " <== Closed file " << fileenv << G4endl;
}
