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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"       

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4FieldManager;
class G4UniformMagField;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    DetectorConstruction();
    ~DetectorConstruction(); 
    G4VPhysicalVolume* Construct() override;
  
    void SetTargetMaterial( const G4String name );
    inline G4Material* GetTargetMaterial() const;
    inline void SetTargetInnerRadius( const G4double value );
    inline G4double GetTargetInnerRadius() const;
    inline void SetTargetOuterRadius( const G4double value );
    inline G4double GetTargetOuterRadius() const;
    void SetMagField( const G4double fieldValue );
    void UpdateGeometry();
  
  private:
    // To be invoked each time the geometry needs to be updated
    G4VPhysicalVolume* ConstructLayer();
  
    void PrintParameters();
    G4Material* fTargetMaterial;
    G4LogicalVolume* fLogicExperimentalHall;
    G4VPhysicalVolume* fPhysExperimentalHall;
    G4LogicalVolume*  fLogicTargetLayer;
    G4VPhysicalVolume* fPhysTargetLayer;
    G4FieldManager* fFieldMgr;
    G4UniformMagField* fUniformMagField; 
    DetectorMessenger* fDetectorMessenger;
    G4double fTargetInnerRadius;
    G4double fTargetOuterRadius;
};

inline G4Material* DetectorConstruction::GetTargetMaterial() const {
  return fTargetMaterial;
}

inline void DetectorConstruction::SetTargetInnerRadius( const G4double value ) {
  fTargetInnerRadius = value;
}

inline G4double DetectorConstruction::GetTargetInnerRadius() const {
  return fTargetInnerRadius;
}

inline void DetectorConstruction::SetTargetOuterRadius( const G4double value ) {
  fTargetOuterRadius = value;
}

inline G4double DetectorConstruction::GetTargetOuterRadius() const {
  return fTargetOuterRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
