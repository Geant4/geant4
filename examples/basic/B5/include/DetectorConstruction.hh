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
/// \file B5/include/DetectorConstruction.hh
/// \brief Definition of the B5::DetectorConstruction class

#ifndef B5DetectorConstruction_h
#define B5DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"

#include <vector>

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;

namespace B5
{

class MagneticField;

/// Detector construction

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    void SetArmAngle(G4double val);
    G4double GetArmAngle() { return fArmAngle; }

    void ConstructMaterials();

  private:
    void DefineCommands();

    G4GenericMessenger* fMessenger;

    static G4ThreadLocal MagneticField* fMagneticField;
    static G4ThreadLocal G4FieldManager* fFieldMgr;

    G4LogicalVolume* fHodoscope1Logical = nullptr;
    G4LogicalVolume* fHodoscope2Logical = nullptr;
    G4LogicalVolume* fWirePlane1Logical = nullptr;
    G4LogicalVolume* fWirePlane2Logical = nullptr;
    G4LogicalVolume* fCellLogical = nullptr;
    G4LogicalVolume* fHadCalScintiLogical = nullptr;
    G4LogicalVolume* fMagneticLogical = nullptr;

    std::vector<G4VisAttributes*> fVisAttributes;

    G4double fArmAngle = 30.*CLHEP::deg;
    G4RotationMatrix* fArmRotation = nullptr;
    G4VPhysicalVolume* fSecondArmPhys = nullptr;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
