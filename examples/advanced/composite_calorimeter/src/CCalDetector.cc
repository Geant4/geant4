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

G4String CCalDetector::pathName = getenv("CCAL_GEOMPATH");

CCalDetector::CCalDetector(const G4String &name):
  detectorName(name) {
#ifdef ddebug
    cout << "CCAL_GEOMPATH=" << pathName << endl;
#endif
    fileName      = 
      CCalGeometryConfiguration::getInstance()->getFileName(name);
    constructFlag = 
      CCalGeometryConfiguration::getInstance()->getConstructFlag(name);
}

CCalDetector::~CCalDetector() {
  CCalDetectorTable::iterator ite;
  for( ite = theDetectorsInside.begin(); ite != theDetectorsInside.end(); ite++ ) {
    delete *ite;
  }
  theDetectorsInside.clear();
}

void CCalDetector::construct() {
#ifdef debug
  cout << "===> Entering CCalDetector::construct() for " << CCalDetectorName << endl;
#endif
  int isgood = 0;

  //If constructFlag is unset we don't go into all this bussines
  if (constructFlag!=0) {
    
    if (!isgood)
      isgood = buildFromFile();
    
    if (isgood) {
      constructDaughters();
      for (unsigned int i=0; i < theDetectorsInside.size(); i++) {
	theDetectorsInside[i]->constructHierarchy();
      }
    }
  }
#ifdef debug
  cout << "===> Exiting CCalDetector::construct() for " << detectorName 
       << endl;
#endif
}


void CCalDetector::addDetector(CCalDetector* det) {
  theDetectorsInside.push_back(det);
}



//========================================================================
//Protected and private methods.
int CCalDetector::buildFromFile() {
  return readFile();
}



//========================================================================
//Global operators
ostream& operator<<(ostream& os, const CCalDetector& det) {
  os << "Detector \"" << det.detectorName 
     << "\" read from " << det.fileName << "." << endl;

  os << "With " << det.theDetectorsInside.size() 
     << " detectors inside { "<< endl;

  for (unsigned int i=0; i<det.theDetectorsInside.size(); i++)
    os << det.theDetectorsInside[i] << endl;

  os << "}" << endl;

  return os;
}







