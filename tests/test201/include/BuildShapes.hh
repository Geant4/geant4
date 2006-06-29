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
//
// $Id: BuildShapes.hh,v 1.4 2006-06-29 21:46:54 gunter Exp $
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
