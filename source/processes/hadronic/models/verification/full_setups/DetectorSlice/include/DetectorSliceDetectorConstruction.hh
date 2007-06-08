#ifndef DetectorSliceDetectorConstruction_H
#define DetectorSliceDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"       

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorSliceDetectorMessenger;
class DetectorSliceSensitiveEmCalo;
class DetectorSliceSensitiveHadCalo;
class G4VisAttributes;


class DetectorSliceDetectorConstruction : public G4VUserDetectorConstruction {

public:

  DetectorSliceDetectorConstruction();
  ~DetectorSliceDetectorConstruction();
  
  G4VPhysicalVolume* Construct();

  void SetTrackerMaterial( const G4String name );
  void SetEmAbsorberMaterial( const G4String name );
  void SetEmActiveMaterial( const G4String name );
  void SetHadAbsorberMaterial( const G4String name );
  void SetHadActiveMaterial( const G4String name );
  void SetMuonMaterial( const G4String name );
  // Use by the messenger.

  inline G4Material* GetTrackerMaterial() const;
  inline G4Material* GetEmAbsorberMaterial() const;
  inline G4Material* GetEmActiveMaterial() const;
  inline G4Material* GetHadAbsorberMaterial() const;
  inline G4Material* GetHadActiveMaterial() const;
  inline G4Material* GetMuonMaterial() const;

  inline void SetTrackerLength( const G4double value );

  inline void SetIsEmCalHomogeneous( const G4bool choice );
  inline void SetEmAbsorberTotalLength( const G4double value );
  inline void SetEmActiveLayerNumber( const G4int value );
  inline void SetEmActiveLayerSize( const G4double value );

  inline void SetIsHadCalHomogeneous( const G4bool choice );
  inline void SetHadAbsorberTotalLength( const G4double value );
  inline void SetHadActiveLayerNumber( const G4int value );
  inline void SetHadActiveLayerSize( const G4double value );

  inline void SetMuonLength( const G4double value );

  inline void SetDetectorRadius( const G4double value );

  void UpdateGeometry();

private:

  void DefineMaterials();
  // Define all the materials.

  G4VPhysicalVolume* ConstructCalorimeter();     
  // To be invoked each time the geometry needs to be updated.

  G4bool areParametersOK();
  // Return true if all the parameters are sensible, false otherwise.

  void PrintParameters();
  // Print the various parameters which define the calorimeter.

  G4Material* Vacuum;
  G4Material* Iron;
  G4Material* Copper;
  G4Material* Tungsten;
  G4Material* Lead;
  G4Material* Uranium;
  G4Material* PbWO4;
  G4Material* Polystyrene;
  G4Material* LiquidArgon;
  G4Material* Silicon;
  G4Material* Quartz;

  G4Material* theTrackerMaterial;
  G4Material* theEmAbsorberMaterial;
  G4Material* theEmActiveMaterial;
  G4Material* theHadAbsorberMaterial;
  G4Material* theHadActiveMaterial;
  G4Material* theMuonMaterial;
  
  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;
  // World envelope. 
  
  G4LogicalVolume*  logicTracker;
  G4VPhysicalVolume* physiTracker;
  // Tracker.

  G4LogicalVolume*  logicEmCalo;
  G4VPhysicalVolume* physiEmCalo;
  // "EM Calorimeter".

  G4LogicalVolume*  logicEmModule;
  G4VPhysicalVolume* physiEmModule;
  // Module of the "EM calorimeter".
  
  G4LogicalVolume*  logicEmAbsorber;
  G4VPhysicalVolume* physiEmAbsorber;
  // Absorber layer of the "EM calorimeter".
  
  G4LogicalVolume*  logicEmActive;
  G4VPhysicalVolume* physiEmActive;
  // Active layer of the "EM calorimeter".

  G4LogicalVolume*  logicHadCalo;
  G4VPhysicalVolume* physiHadCalo;
  // "HAD Calorimeter".

  G4LogicalVolume*  logicHadModule;
  G4VPhysicalVolume* physiHadModule;
  // Module of the "HAD calorimeter".
  
  G4LogicalVolume*  logicHadAbsorber;
  G4VPhysicalVolume* physiHadAbsorber;
  // Absorber layer of the "HAD calorimeter".
  
  G4LogicalVolume*  logicHadActive;
  G4VPhysicalVolume* physiHadActive;
  // Active layer of the "HAD calorimeter".

  G4LogicalVolume*  logicMuon;
  G4VPhysicalVolume* physiMuon;
  // Muon detector.

  DetectorSliceDetectorMessenger* detectorMessenger;
  // Pointer to the Messenger.

  DetectorSliceSensitiveEmCalo* theSensitiveEmCalorimeter;
  DetectorSliceSensitiveHadCalo* theSensitiveHadCalorimeter;
  // Pointers to the sensitive parts of EM and HAD calorimeters..

  G4VisAttributes* theVisAttAbsorber;
  G4VisAttributes* theVisAttActive;
  // Pointers to visual attributes.

  G4bool theIsEmCalHomogeneous; 
  G4bool theIsHadCalHomogeneous; 
  // If "false" then Sampling calorimeter;
  // if "true"  then Homogeneous calorimeter.

  G4double theTrackerLength;

  G4double theEmAbsorberTotalLength;
  G4double theHadAbsorberTotalLength;
  // This is the total length of the absorber material, in [mm],
  // for the electromagnetic and hadronic calorimeters, respectively. 
  // Notice that in the case of a sampling calorimeter, the 
  // active layers are not counted; in the case of an homogenous
  // calorimeter, this length accounts for the overall dimension
  // of the calorimeter.
  // It is assumed that a zero length means that the calorimeter 
  // is not present.

  G4double theMuonLength;

  G4double theDetectorRadius;
  // This is the radius, in [mm], of the detector slice, assumed
  // to be a cylinder.

  G4int theEmActiveLayerNumber;
  G4int theHadActiveLayerNumber;
  G4double theEmActiveLayerSize;
  G4double theHadActiveLayerSize;
  // Number of active layers and length of each of them, in [mm].
  // In the case of a homogeneous calorimeter, these values are
  // neglected.

};


inline G4Material* DetectorSliceDetectorConstruction::
GetTrackerMaterial() const {
  return theTrackerMaterial;
}

inline G4Material* DetectorSliceDetectorConstruction::
GetEmAbsorberMaterial() const {
  return theEmAbsorberMaterial;
}

inline G4Material* DetectorSliceDetectorConstruction::
GetEmActiveMaterial() const {
  return theEmActiveMaterial;
}

inline G4Material* DetectorSliceDetectorConstruction::
GetHadAbsorberMaterial() const {
  return theHadAbsorberMaterial;
}

inline G4Material* DetectorSliceDetectorConstruction::
GetHadActiveMaterial() const {
  return theHadActiveMaterial;
}

inline G4Material* DetectorSliceDetectorConstruction::
GetMuonMaterial() const {
  return theMuonMaterial;
}


inline void DetectorSliceDetectorConstruction::
SetTrackerLength( const G4double value ) {
  theTrackerLength = value;
}

inline void DetectorSliceDetectorConstruction::
SetIsEmCalHomogeneous( const G4bool choice ) {
  theIsEmCalHomogeneous = choice;
}

inline void DetectorSliceDetectorConstruction::
SetEmAbsorberTotalLength( const G4double value ) {
  theEmAbsorberTotalLength = value;
}

inline void DetectorSliceDetectorConstruction::
SetEmActiveLayerNumber( const G4int value ) {
  theEmActiveLayerNumber = value;
}

inline void DetectorSliceDetectorConstruction::
SetEmActiveLayerSize( const G4double value ) {
  theEmActiveLayerSize = value;
}

inline void DetectorSliceDetectorConstruction::
SetIsHadCalHomogeneous( const G4bool choice ) {
  theIsHadCalHomogeneous = choice;
}

inline void DetectorSliceDetectorConstruction::
SetHadAbsorberTotalLength( const G4double value ) {
  theHadAbsorberTotalLength = value;
}

inline void DetectorSliceDetectorConstruction::
SetHadActiveLayerNumber( const G4int value ) {
  theHadActiveLayerNumber = value;
}

inline void DetectorSliceDetectorConstruction::
SetHadActiveLayerSize( const G4double value ) {
  theHadActiveLayerSize = value;
}

inline void DetectorSliceDetectorConstruction::
SetMuonLength( const G4double value ) {
  theMuonLength = value;
}

inline void DetectorSliceDetectorConstruction::
SetDetectorRadius( const G4double value ) {
  theDetectorRadius = value;
}


#endif

