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
// File: CCalDetector.hh
// Description: CCalDetector will provide the basic functionallity of a 
//              detector factory. It has a construct method that describes 
//              the construction sequence of the detector. It has a Name, 
//              a file associated with it, and a collection of other 
//              CCalDetectors that are supposed to be inside this one.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalDetector_h
#define CCalDetector_h 1

#include <iostream>
#include <vector>
#include "globals.hh"

//Forward declartion for the CCalDetectorTable typedef
class CCalDetector;

//A table to hold a list of pointers to CMS Detectors
typedef std::vector<CCalDetector*> CCalDetectorTable;

////////////////////
//At last the class
class CCalDetector {

  friend std::ostream& operator<<(std::ostream&, const CCalDetector&);

public:
  ////////////////////////////////////////////////////////////////
  //Constructor with a name and a filename.
  CCalDetector(const G4String &name);
  ////////////////////////////////////////////////////////////////
  //Destructor
  virtual ~CCalDetector();


  ////////////////////////////////////////////////////////////////
  // Construction related methods

  //This starts the detector construction
  void constructHierarchy() { construct();}
  void construct();


  ////////////////////////////////////////////////////////////////
  //A method that allows to add a new detector inside this one.
  void addDetector(CCalDetector*);
  

  ////////////////////////////////////////////////////////////////
  //Other get methods
  G4String Name() const {return detectorName;}
  G4String baseFileName() const {return fileName;}
  G4String File() const {return fileName+".geom";}
  CCalDetector* getDaughter(G4int i) const {return theDetectorsInside[i];}
  G4int getNDaughters() const {return theDetectorsInside.size();}
  

  ////////////////////////////////////////////////////////////////
  //Local operators
  //Equality only checks name !!!
  G4bool operator==(const CCalDetector& left) const {
    return (detectorName==left.detectorName);
  }
  G4bool operator!=(const CCalDetector& left) const {
    return (detectorName!=left.detectorName);
  }


protected:

  ////////////////////////////////////////////////////////////////
  //Pure Virtual methods  

  //Should read a file and store the information inside the concrete 
  //CCalDetector
  virtual G4int readFile() = 0;
  //Construct the daughters by calling the apropiate constructors
  virtual void constructDaughters() = 0;

  ////////////////////////////////////////////////////////////////
  //Building related methods  

  //Builds the detector from a file
  G4int buildFromFile();

protected:
  ////////////////////////////////////////////////////////////////
  //Data Members

  G4String detectorName;                //Detector name
  G4String fileName;                    //File name from it will be read
  G4String pathName;             //Path in which to look for files

  CCalDetectorTable theDetectorsInside; //A collection of CCalDetectors inside

  G4int constructFlag;                  //True if this detector is to be built
};

#endif



