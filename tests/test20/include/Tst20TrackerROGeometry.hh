// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20TrackerROGeometry.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20TrackerROGeometry  ------
//
// ************************************************************

#ifndef Tst20TrackerROGeometry_h
#define Tst20TrackerROGeometry_h 1

#include "G4VReadOutGeometry.hh"

class Tst20DetectorConstruction;

class Tst20TrackerROGeometry : public G4VReadOutGeometry
{
public:
  Tst20TrackerROGeometry();
  Tst20TrackerROGeometry(G4String);
  Tst20TrackerROGeometry(G4String, Tst20DetectorConstruction*);
  ~Tst20TrackerROGeometry();

private:
  G4VPhysicalVolume* Build();
  Tst20DetectorConstruction* Tst20Detector;
  //pointer to the geometry
};

#endif




