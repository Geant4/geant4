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

    void SetMaterialTracker( const G4String name );
    inline G4Material* GetMaterialTracker() const;
    void SetMaterialEmCalo( const G4String name );
    inline G4Material* GetMaterialEmCalo() const;
    void SetMaterialHadCalo( const G4String name );
    inline G4Material* GetMaterialHadCalo() const;
  
    inline void SetInnerRadiusTracker( const G4double value );
    inline G4double GetInnerRadiusTracker() const;
    inline void SetOuterRadiusTracker( const G4double value );
    inline G4double GetOuterRadiusTracker() const;
    inline void SetInnerRadiusEmCalo( const G4double value );
    inline G4double GetInnerRadiusEmCalo() const;
    inline void SetOuterRadiusEmCalo( const G4double value );
    inline G4double GetOuterRadiusEmCalo() const;
    inline void SetInnerRadiusHadCalo( const G4double value );
    inline G4double GetInnerRadiusHadCalo() const;
    inline void SetOuterRadiusHadCalo( const G4double value );
    inline G4double GetOuterRadiusHadCalo() const;
    inline G4double GetScoringThickness() const;
  
    void UpdateGeometry();
  
  private:
    G4VPhysicalVolume* ConstructDetector();  // To be invoked each time the geometry needs
                                             // to be updated.
    void PrintParameters();
    G4Material* fMaterialTracker;
    G4Material* fMaterialEmCalo;
    G4Material* fMaterialHadCalo;
    G4LogicalVolume*   fExperimentalHall_log;
    G4VPhysicalVolume* fExperimentalHall_phys;
    G4LogicalVolume*   fLogicTrackerShell;
    G4VPhysicalVolume* fPhysiTrackerShell;
    G4LogicalVolume*   fLogicEmCaloShell;
    G4VPhysicalVolume* fPhysiEmCaloShell;
    G4LogicalVolume*   fLogicHadCaloShell;
    G4VPhysicalVolume* fPhysiHadCaloShell;
    G4LogicalVolume*   fLogicScoringTrackerShell;
    G4VPhysicalVolume* fPhysiScoringTrackerShell;
    G4LogicalVolume*   fLogicScoringEmCaloShell;
    G4VPhysicalVolume* fPhysiScoringEmCaloShell;
    G4LogicalVolume*   fLogicScoringHadCaloShell;
    G4VPhysicalVolume* fPhysiScoringHadCaloShell;
    DetectorMessenger* fDetectorMessenger;  
    G4double fInnerRadiusTracker;
    G4double fOuterRadiusTracker;
    G4double fInnerRadiusEmCalo;
    G4double fOuterRadiusEmCalo;
    G4double fInnerRadiusHadCalo;
    G4double fOuterRadiusHadCalo;
    const G4double fScoringThickness = 10.0;  //***LOOKHERE*** thickness of the scoring shell
};

inline G4Material* DetectorConstruction::GetMaterialTracker() const {
  return fMaterialTracker;
}

inline G4Material* DetectorConstruction::GetMaterialEmCalo() const {
  return fMaterialEmCalo;
}

inline G4Material* DetectorConstruction::GetMaterialHadCalo() const {
  return fMaterialHadCalo;
}

inline G4double DetectorConstruction::GetInnerRadiusTracker() const {
  return fInnerRadiusTracker;
}

inline void DetectorConstruction::SetInnerRadiusTracker( const G4double value ) {
  fInnerRadiusTracker = value;
}

inline G4double DetectorConstruction::GetOuterRadiusTracker() const {
  return fOuterRadiusTracker;
}

inline void DetectorConstruction::SetOuterRadiusTracker( const G4double value ) {
  fOuterRadiusTracker = value;
}

inline G4double DetectorConstruction::GetInnerRadiusEmCalo() const {
  return fInnerRadiusEmCalo;
}

inline void DetectorConstruction::SetInnerRadiusEmCalo( const G4double value ) {
  fInnerRadiusEmCalo = value;
}

inline G4double DetectorConstruction::GetOuterRadiusEmCalo() const {
  return fOuterRadiusEmCalo;
}

inline void DetectorConstruction::SetOuterRadiusEmCalo( const G4double value ) {
  fOuterRadiusEmCalo = value;
}

inline G4double DetectorConstruction::GetInnerRadiusHadCalo() const {
  return fInnerRadiusHadCalo;
}

inline void DetectorConstruction::SetInnerRadiusHadCalo( const G4double value ) {
  fInnerRadiusHadCalo = value;
}

inline G4double DetectorConstruction::GetOuterRadiusHadCalo() const {
  return fOuterRadiusHadCalo;
}

inline void DetectorConstruction::SetOuterRadiusHadCalo( const G4double value ) {
  fOuterRadiusHadCalo = value;
}

inline G4double DetectorConstruction::GetScoringThickness() const {
  return fScoringThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
