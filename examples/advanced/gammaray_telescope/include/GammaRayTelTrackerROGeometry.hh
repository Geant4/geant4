// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelTrackerROGeometry.hh,v 1.1 2001-03-05 13:58:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelTrackerROGeometry  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelTrackerROGeometry_h
#define GammaRayTelTrackerROGeometry_h 1

#include "G4VReadOutGeometry.hh"

class GammaRayTelDetectorConstruction;

class GammaRayTelTrackerROGeometry : public G4VReadOutGeometry
{
public:
  GammaRayTelTrackerROGeometry();
  GammaRayTelTrackerROGeometry(G4String);
  GammaRayTelTrackerROGeometry(G4String, GammaRayTelDetectorConstruction*);
  ~GammaRayTelTrackerROGeometry();

private:
  G4VPhysicalVolume* Build();
  GammaRayTelDetectorConstruction* GammaRayTelDetector;
  //pointer to the geometry
};

#endif


