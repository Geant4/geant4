//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: BuildShapes.hh,v 1.3 2001-07-11 10:09:24 gunter Exp $
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
