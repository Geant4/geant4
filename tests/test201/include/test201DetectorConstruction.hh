// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test201DetectorConstruction.hh,v 1.2 1999-12-15 14:55:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Detector Construction for visualization testing.
// John Allison 24th April 1997

#ifndef test201DetectorConstruction_hh
#define test201DetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"

class test201DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  test201DetectorConstruction();
  ~test201DetectorConstruction();

public:
  G4VPhysicalVolume* Construct();
  void SetDetector (G4VPhysicalVolume* pDetector) {fpDetector = pDetector;}

private:
  G4VPhysicalVolume* fpDetector;
};

#endif
