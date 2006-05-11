// $Id: MyDetectorConstruction.hh,v 1.1 2006-05-11 04:35:32 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   MyDetectorConstruction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef MY_DETECTOR_CONSTRUCTION_H
#define MY_DETECTOR_CONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class G4LogicalVolume;
class G4VSensitiveDetector;

class MyDetectorConstruction : public G4VUserDetectorConstruction {
private:
  G4LogicalVolume* scoreVoxel;

public:
  MyDetectorConstruction();
  ~MyDetectorConstruction(); 

  virtual G4VPhysicalVolume* Construct();
  
  void SetSDtoScoreVoxel(G4VSensitiveDetector* asd);

};

#endif
