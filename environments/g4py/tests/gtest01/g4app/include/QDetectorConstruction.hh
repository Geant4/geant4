// $Id: QDetectorConstruction.hh,v 1.1 2006-02-27 10:05:24 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   QDetectorConstruction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef Q_DETECTOR_CONSTRUCTION_H
#define Q_DETECTOR_CONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class QDetectorConstruction : public G4VUserDetectorConstruction {

public:
  QDetectorConstruction();
  ~QDetectorConstruction(); 

  virtual G4VPhysicalVolume* Construct(); 

};

#endif
