///////////////////////////////////////////////////////////////////////////////
// File: CMSDetector.hh
// Date: 03/98 
// Last modification: 24/03/98 I.G.
//                    03/05/99 I.G. Added sensitive behaviour
//                    04/05/99 I.G. Move to use GeometryConfiguration.
//                    16/12/99 I.G. RWTPtrOrderedVector ->G4RWTPtrOrderedVector
//                    17/11/01 P.A. G4RWTPtrOrderedVector ->STL vector
//                                  Needs change to STL!!!
//                    15/05/02 S.B. Clean up persistency code
// Description: CMSDetector will provide the basic functionallity of a CMS
//              "detector". Based on the ideas from Sunanda and myself, it
//              has a construct method that describes the construction
//              sequence of the detector. It has a Name, a file associated
//              with it (objectivity can also be used), and
//              a collection of other CMSDetectors that are supposed
//              to be inside this one.
///////////////////////////////////////////////////////////////////////////////
#ifndef CMSDetector_h
#define CMSDetector_h 1

#include <iostream>
#include <vector>
#include "globals.hh"

//Forward declartion for the CMSDetectorTable typedef
class CMSDetector;

//A table to hold a list of pointers to CMS Detectors
typedef  vector<CMSDetector*> CMSDetectorTable;

////////////////////
//At last the class
class CMSDetector {

  friend ostream& operator<<(ostream&, const CMSDetector&);

public:
  ////////////////////////////////////////////////////////////////
  //Constructor with a name and a filename.
  CMSDetector(const G4String &name);
  ////////////////////////////////////////////////////////////////
  //Destructor
  virtual ~CMSDetector();


  ////////////////////////////////////////////////////////////////
  // Construction related methods

  //This starts the detector construction
  void constructHierarchy() { construct();}
  void construct();


  ////////////////////////////////////////////////////////////////
  //A method that allows to add a new detector inside this one.
  void addDetector(CMSDetector*);
  

  ////////////////////////////////////////////////////////////////
  //Other get methods
  G4String Name() const {return cmsDetectorName;}
  G4String baseFileName() const {return fileName;}
  G4String File() const {return fileName+".geom";}
  CMSDetector* getDaughter(int i) const {return theDetectorsInside[i];}
  int getNDaughters() const {return theDetectorsInside.size();}
  

  ////////////////////////////////////////////////////////////////
  //Local operators
  //Equality only checks name !!!
  G4bool operator==(const CMSDetector& left) const {
    return (cmsDetectorName==left.cmsDetectorName);
  }
  G4bool operator!=(const CMSDetector& left) const {
    return (cmsDetectorName!=left.cmsDetectorName);
  }





protected:

  ////////////////////////////////////////////////////////////////
  //Pure Virtual methods  

  //Should read a file and store the information inside the concrete 
  //CMSDetector
  virtual int readFile() = 0;
  //Construct the daughters by calling the apropiate constructors
  virtual void constructDaughters() = 0;





  ////////////////////////////////////////////////////////////////
  //Building related methods  

  //Builds the detector from a file
  int buildFromFile();

#ifdef CMS_USE_OBJY
  //Builds the detector from objectivity
  virtual int buildFromPersistent() {
    cout << "WARNING: Persistent to Transient not yet implemented for " 
	 << cmsDetectorName << endl;
    return 0;
  }

  //Make this detector persistent
  virtual ooHandle(PCMSDetector) makePersistent() const {
    cout << "WARNING: Transient to Persistent not yet implemented for "
	 << cmsDetectorName << endl;
    return 0;
  }
#endif




protected:

  ////////////////////////////////////////////////////////////////
  //Data Members

  G4String cmsDetectorName;             //Detector name
  G4String fileName;                    //File name from it will be read
  static G4String pathName;             //Path in which to look for files
  
  CMSDetectorTable theDetectorsInside;  //A collection of CMSDetectors inside

  int constructFlag;                    //True if this detector is to be built
};

#endif



