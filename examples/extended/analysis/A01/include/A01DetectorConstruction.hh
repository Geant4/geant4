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
/// \file analysis/A01/include/A01DetectorConstruction.hh
/// \brief Definition of the A01DetectorConstruction class
//
// $Id$
// --------------------------------------------------------------
//

#ifndef A01DetectorConstruction_h
#define A01DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"

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

    virtual G4VPhysicalVolume* Construct();
    void SetArmAngle(G4double val);
    inline G4double GetArmAngle() { return fArmAngle; }

    void ConstructMaterials();
    void DestroyMaterials();
    void DumpGeometricalTree(G4VPhysicalVolume* aVolume,G4int depth=0);

  private:
    A01DetectorConstMessenger* fMessenger;
    A01MagneticField* fMagneticField;
    G4FieldManager* fFieldMgr;

    G4Material* fAir;
    G4Material* fArgonGas;
    G4Material* fScintillator;
    G4Material* fCsI;
    G4Material* fLead;

    G4VisAttributes* fWorldVisAtt;
    G4VisAttributes* fMagneticVisAtt;
    G4VisAttributes* fArmVisAtt;
    G4VisAttributes* fHodoscopeVisAtt;
    G4VisAttributes* fChamberVisAtt;
    G4VisAttributes* fWirePlaneVisAtt;
    G4VisAttributes* fEMcalorimeterVisAtt;
    G4VisAttributes* fCellVisAtt;
    G4VisAttributes* fHadCalorimeterVisAtt;
    G4VisAttributes* fHadCalorimeterCellVisAtt;

    G4double fArmAngle;
    G4RotationMatrix* fArmRotation;
    G4VPhysicalVolume* fSecondArmPhys;
};

#endif

