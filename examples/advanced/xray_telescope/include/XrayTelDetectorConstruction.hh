// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelDetectorConstruction.hh                  *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of X-ray Telescope advanced example.
// - Based on Chandra and XMM models
//
//
// **********************************************************************

#ifndef XrayTelDetectorConstruction_H
#define XrayTelDetectorConstruction_H 1

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"

class XrayTelDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  XrayTelDetectorConstruction();
  ~XrayTelDetectorConstruction();

public:
  G4VPhysicalVolume* Construct();

private:
  G4double world_x;
  G4double world_y;
  G4double world_z;
  void ConstructTelescope();
  void ConstructFocalPlane();
  G4VPhysicalVolume* physicalWorld;
};

#endif

