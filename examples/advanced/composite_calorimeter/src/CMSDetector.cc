///////////////////////////////////////////////////////////////////////////////
// File: CMSDetector.cc
// Date: 03/98 I. Gonzalez
// Modifications: 24/03/98 I.G.
//                03/05/99 I.G. Added sensitive behaviour
//                02/05/00 I.G. S.B. Use local file if available
//                15/05/02 S.B. Remove persistency
///////////////////////////////////////////////////////////////////////////////

#include "CMSDetector.hh"

#include <fstream>
#include "GeometryConfiguration.hh"
#include "utils.hh"

//Comment/Uncomment to hide/show some debug information
//#define debug
//Comment/Uncomment to hide/show some more debug information
//#define ddebug

G4String CMSDetector::pathName = getenv("OSCARGEOMPATH");

CMSDetector::CMSDetector(const G4String &name):
  cmsDetectorName(name) {
#ifdef ddebug
    cout << "OSCARGEOMPATH=" << pathName << endl;
#endif
    fileName      = 
      GeometryConfiguration::getInstance()->getFileName(name);
    constructFlag = 
      GeometryConfiguration::getInstance()->getConstructFlag(name);
}

CMSDetector::~CMSDetector() {
  CMSDetectorTable::iterator ite;
  for( ite = theDetectorsInside.begin(); ite != theDetectorsInside.end(); ite++ ) {
    delete *ite;
  }
  theDetectorsInside.clear();
}

void CMSDetector::construct() {
#ifdef debug
  cout << "===> Entering CMSDetector::construct() for " << cmsDetectorName << endl;
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
  cout << "===> Exiting CMSDetector::construct() for " << cmsDetectorName 
       << endl;
#endif
}


void CMSDetector::addDetector(CMSDetector* det) {
  theDetectorsInside.push_back(det);
}



//========================================================================
//Protected and private methods.
int CMSDetector::buildFromFile() {
  return readFile();
}



//========================================================================
//Global operators
ostream& operator<<(ostream& os, const CMSDetector& det) {
  os << "Detector \"" << det.cmsDetectorName 
     << "\" read from " << det.fileName << "." << endl;

  os << "With " << det.theDetectorsInside.size() 
     << " detectors inside { "<< endl;

  for (unsigned int i=0; i<det.theDetectorsInside.size(); i++)
    os << det.theDetectorsInside[i] << endl;

  os << "}" << endl;

  return os;
}







