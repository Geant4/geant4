// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelPayloadROGeometry.hh,v 1.1 2000-11-15 20:27:39 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelPayloadROGeometry  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelPayloadROGeometry_h
#define GammaRayTelPayloadROGeometry_h 1

#include "G4VReadOutGeometry.hh"

class GammaRayTelDetectorConstruction;

class GammaRayTelPayloadROGeometry : public G4VReadOutGeometry
{
public:
  GammaRayTelPayloadROGeometry();
  GammaRayTelPayloadROGeometry(G4String);
  GammaRayTelPayloadROGeometry(G4String, GammaRayTelDetectorConstruction*);
  ~GammaRayTelPayloadROGeometry();

private:
  G4VPhysicalVolume* Build();
  GammaRayTelDetectorConstruction* GammaRayTelDetector;
  //pointer to the geometry
};

#endif


