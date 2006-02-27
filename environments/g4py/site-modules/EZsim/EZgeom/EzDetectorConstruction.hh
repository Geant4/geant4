// $Id: EzDetectorConstruction.hh,v 1.1 2006-02-27 09:46:31 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   EzDetectorConstruction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef EZ_DETECTOR_CONSTRUCTION_H
#define EZ_DETECTOR_CONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class EzDetectorConstruction : public G4VUserDetectorConstruction {

public:
  EzDetectorConstruction();
  ~EzDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

};

#endif
