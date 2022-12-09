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
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    DetectorConstruction();
    ~DetectorConstruction();
    G4VPhysicalVolume* Construct();
    void SetMaterial( const G4String name );
    inline G4Material* GetMaterial() const;
    inline void SetRadius( const G4double value );
    inline G4double GetRadius() const;
    void UpdateGeometry();
  private:
    G4VPhysicalVolume* ConstructSphere();  // To be invoked each time the geometry needs
                                           // to be updated
    void PrintParameters();
    G4Material* fMaterial;
    G4LogicalVolume*   fExperimentalHall_log;
    G4VPhysicalVolume* fExperimentalHall_phys;
    G4LogicalVolume*   fLogicSphere;
    G4VPhysicalVolume* fPhysiSphere;
    G4LogicalVolume*   fLogicScoringShell;
    G4VPhysicalVolume* fPhysiScoringShell;
    DetectorMessenger* fDetectorMessenger;
    G4double fRadius;
    const G4double fScoringThickness = 10.0;  //***LOOKHERE*** thickness of the scoring shell
};

inline G4Material* DetectorConstruction::GetMaterial() const {
  return fMaterial;
}

inline void DetectorConstruction::SetRadius( const G4double value ) {
  fRadius = value;
}

inline G4double DetectorConstruction::GetRadius() const {
  return fRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
