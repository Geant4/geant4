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
#ifndef Tst18DetectorConstruction_h
#define Tst18DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4VPhysicalVolume.hh"   

class Tst18DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst18DetectorConstruction();
    ~Tst18DetectorConstruction();
    void AddBox (G4String prefix,
      G4ThreeVector p1, G4ThreeVector p2, G4ThreeVector p3,
      G4ThreeVector p4, G4ThreeVector p5, G4ThreeVector p6,
      G4ThreeVector p7, G4ThreeVector p8,
      G4Material* theMaterial, G4VPhysicalVolume* theParentVolume);
    void AddGeometry001 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry002 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry003 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry004 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry005 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry006 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry007 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry008 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry009 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry010 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry011 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry012 (G4Material* Ni, G4Material* Au, G4Material* Invar, G4VPhysicalVolume* theParentVolume);
    void AddGeometry021 (G4Material* Cf, G4Material* Epoxy, G4Material* Au, G4VPhysicalVolume* theParentVolume);
    void AddGeometry022 (G4Material* Cf, G4Material* Epoxy, G4Material* Au, G4VPhysicalVolume* theParentVolume);
    void AddGeometry023 (G4Material* Cf, G4Material* Epoxy, G4Material* Au, G4VPhysicalVolume* theParentVolume);
    void AddGeometry024 (G4Material* Cf, G4Material* Epoxy, G4Material* Au, G4VPhysicalVolume* theParentVolume);
    void AddGeometry025 (G4Material* Cf, G4Material* Epoxy, G4Material* Au, G4VPhysicalVolume* theParentVolume);
    void AddGeometry026 (G4Material* Cf, G4Material* Epoxy, G4Material* Au, G4VPhysicalVolume* theParentVolume);
    void AddGeometry027 (G4Material* Cf, G4Material* Epoxy, G4Material* Au, G4VPhysicalVolume* theParentVolume);
    void AddGeometry028 (G4Material* Cf, G4Material* Epoxy, G4Material* Au, G4VPhysicalVolume* theParentVolume);
    void AddGeometry029 (G4Material* Cf, G4Material* Epoxy, G4Material* Au, G4VPhysicalVolume* theParentVolume);


  public:
     G4VPhysicalVolume* Construct();

  private:
     G4double Spacecraft_x;
     G4double Spacecraft_y;
     G4double Spacecraft_z;

};

#endif

