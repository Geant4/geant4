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
#ifndef DicomBeamDevice__HH
#define DicomBeamDevice__HH

#include "dcmtk/dcmdata/dcfilefo.h"

#include <vector>
#include <iostream>
class DicomBeamDeviceDevice;
class DicomBeamDeviceControlPoint;

#include "dcmtk/dcmrt/seq/drtblds1.h" // for BeamLimitingDeviceSequenceInRTBeamsModule
#include "dcmtk/dcmrt/seq/drtbldps.h"  // for BeamLimitingDevicePositionSequence

class DicomBeamDevice 
{ 
public:
  DicomBeamDevice(DRTBeamLimitingDeviceSequenceInRTBeamsModule::Item bldsItem);
  DicomBeamDevice(DRTBeamLimitingDevicePositionSequence::Item bldpsItem);
  ~DicomBeamDevice(){};

public:
  void SetSourceToBeamLimitingDeviceDistance(Float64 dat){
    theSourceToBeamLimitingDeviceDistance= dat;
  }
  Float64 GetSourceToBeamLimitingDeviceDistance() const {
    return theSourceToBeamLimitingDeviceDistance;
  }
  void SetNumberOfLeafJawPairs(Sint32 dat){ 
    theNumberOfLeafJawPairs= dat;
  }
  Sint32 GetNumberOfLeafJawPairs() const {
    return theNumberOfLeafJawPairs;
  }
  void SetType(OFString dat){ 
    theType = dat;
  }
  OFString GetType() const {
    return theType;
  }
  void AddPositionBoundary( Float64 dat ){
    thePositionBoundaries.push_back(dat);
  }
  Float64 GetPositionBoundary( size_t ii ) {
    return thePositionBoundaries[ii];
  }
  void Print( std::ostream& out );

private:
  OFString theType;
  Float64 theSourceToBeamLimitingDeviceDistance;
  Sint32  theNumberOfLeafJawPairs;
  std::vector<Float64> thePositionBoundaries;
};

#endif
