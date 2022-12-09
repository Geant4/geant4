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
class G4FieldManager;
class G4UniformMagField;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction {
  public:  
    DetectorConstruction();
    ~DetectorConstruction();
    
    G4VPhysicalVolume* Construct();
    void ConstructSDandField();
  
    void SetMagField( const G4double fieldValue );
    void SetAbsorberMaterial( const G4String name );
    void SetActiveMaterial( const G4String name );
    // Use by the messenger.
  
    inline G4Material* GetAbsorberMaterial() const;
    inline G4Material* GetActiveMaterial() const;
  
    inline void SetIsCalHomogeneous( const G4bool choice );
    inline void SetIsUnitInLambda( const G4bool choice );
    inline void SetAbsorberTotalLength( const G4double value );
    inline void SetCalorimeterRadius( const G4double value );
    inline void SetActiveLayerNumber( const G4int value );
    inline void SetActiveLayerSize( const G4double value );
    // To define the calorimeter geometry.
  
    inline void SetIsRadiusUnitInLambda( const G4bool choice );
    
    void UpdateGeometry();

    inline G4double GetCaloLength() const;

  private:
  
    void DefineMaterials();
    // Define all the materials.
  
    G4VPhysicalVolume* ConstructCalorimeter();     
    // To be invoked each time the geometry needs to be updated.
  
    G4bool AreParametersOK();
    // Return true if all the parameters are sensible, false otherwise.
  
    void PrintParameters();
    // Print the various parameters which define the calorimeter.
  
    G4Material* fVacuum;
    G4Material* fIron;
    G4Material* fCopper;
    G4Material* fTungsten;
    G4Material* fLead;
    G4Material* fUranium;
    G4Material* fPbWO4;
    G4Material* fPolystyrene;
    G4Material* fLiquidArgon;
    G4Material* fSilicon;
    G4Material* fQuartz;
    G4Material* fBrass;
    G4Material* fAluminium;
    G4Material* fGraphite;
    G4Material* fAbsorberMaterial;
    G4Material* fActiveMaterial;
    
    G4LogicalVolume* fExperimentalHall_log;
    G4VPhysicalVolume* fExperimentalHall_phys;
    // World envelope. 
    
    G4LogicalVolume*  fLogicCalo;
    G4VPhysicalVolume* fPhysiCalo;
    // "Calorimeter".
    
    G4LogicalVolume*  fLogicModule;
    G4VPhysicalVolume* fPhysiModule;
    // Module of the "calorimeter".
    
    G4LogicalVolume*  fLogicAbsorber;
    G4VPhysicalVolume* fPhysiAbsorber;
    // Absorber layer of the "calorimeter".
    
    G4LogicalVolume*  fLogicActive;
    G4VPhysicalVolume* fPhysiActive;
    // Active layer of the "calorimeter".
  
    G4FieldManager* fFieldMgr;
    // Pointer to the field manager.
  
    G4UniformMagField* fUniformMagField; 
    // Pointer to the uniform magnetic field.
    
    DetectorMessenger* fDetectorMessenger;
    // Pointer to the Messenger.
  
    G4bool fIsCalHomogeneous; 
    // If false then Sampling calorimeter;
    // If true  then Homogeneous calorimeter.
  
    G4bool fIsUnitInLambda;
    // If false then normal unit of length to express the absorber total length.
    // If true  then lambda (interaction length) to express the absorber total length.
  
    G4double fAbsorberTotalLength;
    // This is the total length of the absorber material, expressed 
    // in unit of length (e.g. m, cm, mm) if theIsUnitInLambda is false, 
    // otherwise in number of lambdas (interaction lengths).
    // Notice that in the case of a sampling calorimeter (i.e. 
    // theIsCalHomogeneous is false), the active layers are not counted;
    // in the case of an homogenous calorimeter, this length account
    // for the overall dimension of the calorimeter.
  
    G4double fCalorimeterRadius;
    // This is the radius of the calorimeter which is a cylinder, expressed 
    // in unit of length (e.g. m, cm, mm) if theIsUnitInLambda is false, 
    // otherwise in number of lambdas (interaction lengths) of the absorber.
  
    G4int fActiveLayerNumber;
    G4double fActiveLayerSize;
    // Number of active layers and length of each of them (in normal unit
    // of length, e.g. mm): in the case of sampling calorimeter
    // (i.e. theIsCalHomogeneous is false) the medium is theActiveMaterial;
    // in the case of an homogeneous calorimeter, the "active layers" are
    // only a fictitious way to sample the longitudinal energy deposits,
    // but they are actually made of the same absorber material, and their
    // thickness is taken into account in theAbsorberTotalLength.
  
    G4bool fIsRadiusUnitInLambda;
    // If false then normal unit of length to express the radius bin size.
    // If true  then lambda (interaction length of the absorber) to express 
    // the radius bin size.
  
    G4double fCaloLength;  // total length of the calorimeter along its (z) axis
  
    // Scoring part
    G4LogicalVolume*  fLogicScoringUpDown;
    G4VPhysicalVolume* fPhysiScoringUpstream;
    G4VPhysicalVolume* fPhysiScoringDownstream;
    G4LogicalVolume*  fLogicScoringSide;
    G4VPhysicalVolume* fPhysiScoringSide;
    const G4double fScoringThickness = 10.0;
};

inline G4Material* DetectorConstruction::GetAbsorberMaterial() const {
  return fAbsorberMaterial;
}

inline G4Material* DetectorConstruction::GetActiveMaterial() const {
  return fActiveMaterial;
}

inline void DetectorConstruction::SetIsCalHomogeneous( const G4bool choice ) {
  fIsCalHomogeneous = choice;
}

inline void DetectorConstruction::SetIsUnitInLambda( const G4bool choice ) {
  fIsUnitInLambda = choice;
}

inline void DetectorConstruction::SetAbsorberTotalLength( const G4double value ) {
  fAbsorberTotalLength = value;
}

inline void DetectorConstruction::SetCalorimeterRadius( const G4double value ) {
  fCalorimeterRadius = value;
}

inline void DetectorConstruction::SetActiveLayerNumber( const G4int value ) {
  fActiveLayerNumber = value;
}

inline void DetectorConstruction::SetActiveLayerSize( const G4double value ) {
  fActiveLayerSize = value;
}

inline void DetectorConstruction::SetIsRadiusUnitInLambda( const G4bool choice ) {
  fIsRadiusUnitInLambda = choice;
}

inline G4double DetectorConstruction::GetCaloLength() const {
  return fCaloLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
