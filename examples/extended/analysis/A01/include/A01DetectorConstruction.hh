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
// $Id: A01DetectorConstruction.hh,v 1.3 2002-12-13 11:34:28 gunter Exp $
// --------------------------------------------------------------
//

#ifndef A01DetectorConstruction_h
#define A01DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class A01DetectorConstMessenger;
class A01MagneticField;

class A01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    A01DetectorConstruction();
    virtual ~A01DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    void SetArmAngle(G4double val);
    inline G4double GetArmAngle() { return armAngle; }

  private:
    void ConstructMaterials();
    void DestroyMaterials();
    void DumpGeometricalTree(G4VPhysicalVolume* aVolume,G4int depth=0);

  private:
    A01DetectorConstMessenger* messenger;
    A01MagneticField* magneticField;

    G4Material* air;
    G4Material* argonGas;
    G4Material* scintillator;
    G4Material* CsI;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    G4Material* lead;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    G4VSensitiveDetector* hodoscope1;
    G4VSensitiveDetector* hodoscope2;
    G4VSensitiveDetector* chamber1;
    G4VSensitiveDetector* chamber2;
    G4VSensitiveDetector* EMcalorimeter;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    G4VSensitiveDetector* HadCalorimeter;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    G4VisAttributes* worldVisAtt;
    G4VisAttributes* magneticVisAtt;
    G4VisAttributes* armVisAtt;
    G4VisAttributes* hodoscopeVisAtt;
    G4VisAttributes* chamberVisAtt;
    G4VisAttributes* wirePlaneVisAtt;
    G4VisAttributes* EMcalorimeterVisAtt;
    G4VisAttributes* cellVisAtt;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    G4VisAttributes* HadCalorimeterVisAtt;
    G4VisAttributes* HadCalorimeterCellVisAtt;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    G4double armAngle;
    G4RotationMatrix* armRotation;
    G4VPhysicalVolume* secondArmPhys;
};

#endif

