// $Id: ExN01DetectorConstruction.hh,v 1.1 2006-02-27 09:47:41 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   ExN01DetectorConstruction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef EXN01_DETECTOR_CONSTRUCTION_H
#define EXN01_DETECTOR_CONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class ExN01DetectorConstruction : public G4VUserDetectorConstruction {
public:
  
  ExN01DetectorConstruction();
  ~ExN01DetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  
private:
  // Logical volumes
  //
  G4LogicalVolume* experimentalHall_log;
  G4LogicalVolume* tracker_log;
  G4LogicalVolume* calorimeterBlock_log;
  G4LogicalVolume* calorimeterLayer_log;
  
  // Physical volumes
  //
  G4VPhysicalVolume* experimentalHall_phys;
  G4VPhysicalVolume* calorimeterLayer_phys;
  G4VPhysicalVolume* calorimeterBlock_phys;
  G4VPhysicalVolume* tracker_phys;

};

#endif

