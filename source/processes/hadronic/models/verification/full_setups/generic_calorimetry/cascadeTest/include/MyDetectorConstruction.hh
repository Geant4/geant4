#ifndef MyDetectorConstruction_H
#define MyDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"       

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UniformMagField;
class G4Material;
class MyDetectorMessenger;


class MyDetectorConstruction : public G4VUserDetectorConstruction {

public:

  MyDetectorConstruction();
  ~MyDetectorConstruction();
  
  G4VPhysicalVolume* Construct();

  void SetMagField(G4double fieldValue);
  void SetAbsorberMaterial(const G4String name);
  void SetActiveMaterial(const G4String name);
  // Use by the messenger.

  inline G4Material* GetAbsorberMaterial() const;
  inline G4Material* GetActiveMaterial() const;

private:

  void PrintParameters();
  // Print the name of absorber and active materials.

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
  G4Material* theAbsorberMaterial;
  G4Material* theActiveMaterial;
  
  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;
  // World envelope. 
  
  G4LogicalVolume*  logicCalo;
  G4VPhysicalVolume* physiCalo;
  // "Calorimeter".
  
  G4LogicalVolume*  logicModule;
  G4VPhysicalVolume* physiModule;
  // Module of the "calorimeter".
  
  G4LogicalVolume*  logicAbsorber;
  G4VPhysicalVolume* physiAbsorber;
  // Absorber layer of the "calorimeter".
  
  G4LogicalVolume*  logicActive;
  G4VPhysicalVolume* physiActive;
  // Active layer of the "calorimeter".

  G4UniformMagField* uniformMagField; 
  // Pointer to the uniform magnetic field.
  
  MyDetectorMessenger* detectorMessenger;
  // pointer to the Messenger
  
};


inline G4Material* MyDetectorConstruction::GetAbsorberMaterial() const {
  return theAbsorberMaterial;
}


inline G4Material* MyDetectorConstruction::GetActiveMaterial() const {
  return theActiveMaterial;
}


#endif

