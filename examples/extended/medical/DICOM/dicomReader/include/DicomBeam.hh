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
#ifndef DicomBeam__HH
#define DicomBeam__HH

#include "dcmtk/dcmdata/dcfilefo.h"
#include "G4ThreeVector.hh"

#include <vector>
#include <iostream>
class DicomVBeamDevice;
class DicomBeamControlPoint;
class DicomBeamCompensator;
class DicomBeamBlock;
class DicomBeamWedge;

class DicomBeam 
{ 
public:
  DicomBeam();
  ~DicomBeam(){};

public:
  void SetDoseSpecificationPoint(G4ThreeVector point ){
    theDoseSpecificationPoint = point;
  }
  void SetMeterset(Float64 dat){
    theMeterset = dat;
  }
  void SetSourceAxisDistance(Float64 dat) {
    theSourceAxisDistance = dat;
  }
  void SetNumber(Sint32 dat){ 
    theNumber = dat;
  }
  void SetRadiationType(OFString dat){ 
    theRadiationType = dat;
  }
  void AddDevice( DicomVBeamDevice* db ){
    theDevices.push_back(db);
  }
  void AddControlPoint( DicomBeamControlPoint* db ){
    theControlPoints.push_back(db);
  }
  void AddCompensator( DicomBeamCompensator* db ){
    theCompensators.push_back(db);
  }
  void AddBlock( DicomBeamBlock* db ){
    theBlocks.push_back(db);
  }
  void AddWedge( DicomBeamWedge* db ){
    theWedges.push_back(db);
  }
  size_t GetNControlPoints() const {
    return theControlPoints.size();
  }
  DicomBeamControlPoint* GetControlPoint( size_t ii ) {
    return theControlPoints[ii];
  }

  void SetControlPointMetersets();
  
  void Print( std::ostream& out );

  void DumpToFile();
  
private:
  G4ThreeVector theDoseSpecificationPoint;
  Float64 theMeterset;
  Float64 theSourceAxisDistance;
  Sint32 theNumber;
  OFString theRadiationType;
  std::vector<DicomVBeamDevice*> theDevices;
  std::vector<DicomBeamControlPoint*> theControlPoints;
  std::vector<DicomBeamCompensator*> theCompensators; 
  std::vector<DicomBeamBlock*> theBlocks; 
  std::vector<DicomBeamWedge*> theWedges; 

};

#endif  
