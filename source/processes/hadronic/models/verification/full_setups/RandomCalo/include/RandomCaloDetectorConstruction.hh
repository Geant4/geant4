#ifndef RandomCaloDetectorConstruction_H
#define RandomCaloDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4FieldManager;
class G4UniformMagField;
class G4Material;
class RandomCaloDetectorMessenger;
class G4NistManager;


class RandomCaloDetectorConstruction : public G4VUserDetectorConstruction {

public:

  RandomCaloDetectorConstruction();
  ~RandomCaloDetectorConstruction();
  
  G4VPhysicalVolume* Construct();

  void SetMagField( const G4double fieldValue );
  // Use by the messenger.

  void UpdateGeometry();

private:

  void DefineMaterials();
  // Define all the materials.

  G4VPhysicalVolume* ConstructCalorimeter();     
  // To be invoked each time the geometry needs to be updated.

  void PrintParameters();
  // Print the various parameters which define the calorimeter.

  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;
  // World envelope. 
  
  G4FieldManager* fieldMgr;
  // Pointer to the field manager.

  G4UniformMagField* uniformMagField; 
  // Pointer to the uniform magnetic field.
  
  RandomCaloDetectorMessenger* detectorMessenger;
  // Pointer to the Messenger.

  G4NistManager* nistManager;
  // Pointer to the NIST material manager.

  std::vector< G4Material* > vecNistMaterials;
  // Vector of pointers to the NIST materials that we consider
  // (i.e. all materials whose elements have Z <= 92 and at
  //  least one stable isotope).

  const G4double thicknessLayer;  // Thickness of each layer.
  const G4double radiusCalo;      // Radius of the cylindrical calorimeter

};


#endif

