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
// File: CCalDetector.cc
// Description: CCalDetector is a detector factory class
///////////////////////////////////////////////////////////////////////////////

#include "CCalDetector.hh"

#include <fstream>
#include "CCalGeometryConfiguration.hh"
#include "CCalutils.hh"

//Comment/Uncomment to hide/show some debug information
//#define debug
//Comment/Uncomment to hide/show some more debug information
//#define ddebug

CCalDetector::CCalDetector(const G4String &name): detectorName(name) 
{
  if (std::getenv("CCAL_GEOMPATH"))    
      pathName = std::getenv("CCAL_GEOMPATH");
  else
    G4Exception("CCalDetector::CCalDetector","ccal001",
                FatalException,
                "Environment variable CCAL_GEOMPATH not defined");
      
#ifdef debug
    G4cout << "CCAL_GEOMPATH=" << pathName << G4endl;
#endif
    fileName      = 
      CCalGeometryConfiguration::getInstance()->getFileName(name);
    constructFlag = 
      CCalGeometryConfiguration::getInstance()->getConstructFlag(name);
}

CCalDetector::~CCalDetector() {
  for (auto ite=theDetectorsInside.begin(); ite !=theDetectorsInside.end(); ++ite) {
    delete *ite;
  }
  theDetectorsInside.clear();
}

void CCalDetector::construct() {
#ifdef debug
  G4cout << "===> Entering CCalDetector::construct() for " << Name() << G4endl;
#endif
  G4int isgood = 0;

  //If constructFlag is unset we don't go into all this bussines
  if (constructFlag!=0) {
    
    if (!isgood)
      isgood = buildFromFile();
    
    if (isgood) {
      constructDaughters();
      for (std::size_t i=0; i < theDetectorsInside.size(); ++i) {
        theDetectorsInside[i]->constructHierarchy();
      }
    }
  }
#ifdef debug
  G4cout << "===> Exiting CCalDetector::construct() for " << Name() << G4endl;
#endif
}

void CCalDetector::addDetector(CCalDetector* det) {
  theDetectorsInside.push_back(det);
}

//========================================================================
//Protected and private methods.
G4int CCalDetector::buildFromFile() {
  return readFile();
}

//========================================================================
//Global operators
std::ostream& operator<<(std::ostream& os, const CCalDetector& det) {
  os << "Detector \"" << det.detectorName 
     << "\" read from " << det.fileName << "." << G4endl;

  os << "With " << det.theDetectorsInside.size() 
     << " detectors inside { "<< G4endl;

  for (std::size_t i=0; i<det.theDetectorsInside.size(); ++i)
    os << det.theDetectorsInside[i] << G4endl;

  os << "}" << G4endl;

  return os;
}
