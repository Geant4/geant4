// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: BuildShapes.hh,v 1.1 1999-04-16 10:32:27 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef BUILDSHAPES_HH
#define BUILDSHAPES_HH

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;

G4VPhysicalVolume* BuildBox();
G4VPhysicalVolume* BuildCylinder();
G4VPhysicalVolume* BuildTubs();
G4VPhysicalVolume* BuildCons();
G4VPhysicalVolume* BuildTrd();
G4VPhysicalVolume* BuildTrap();
G4VPhysicalVolume* BuildSphereFull();
G4VPhysicalVolume* BuildSphereSeg();
G4VPhysicalVolume* BuildPara();
G4VPhysicalVolume* BuildPCon();
G4VPhysicalVolume* BuildPGon();
G4VPhysicalVolume* BuildForcedWireframeBox();

class testshapesDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  testshapesDetectorConstruction(): fpDetector (0) {}
  ~testshapesDetectorConstruction() {}

public:
  G4VPhysicalVolume* Construct() {return fpDetector;}
  void SetDetector (G4VPhysicalVolume* pDetector) {fpDetector = pDetector;}

private:
  G4VPhysicalVolume* fpDetector;
};

#endif
