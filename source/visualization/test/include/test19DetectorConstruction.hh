// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test19DetectorConstruction.hh,v 1.1 1999-04-16 10:32:31 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Detector Construction for visualization testing.
// John Allison 24th April 1997

#ifndef test19DetectorConstruction_hh
#define test19DetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"

class test19DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  test19DetectorConstruction();
  ~test19DetectorConstruction();

public:
  G4VPhysicalVolume* Construct();
  void SetDetector (G4VPhysicalVolume* pDetector) {fpDetector = pDetector;}

private:
  G4VPhysicalVolume* fpDetector;
};

#endif
