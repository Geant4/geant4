#ifndef MyDetectorConstruction_H
#define MyDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"

#include "G4GDMLParser.hh"

class G4VPhysicalVolume;
class G4FieldManager;
class G4UniformMagField;
class G4Material;
class MyDetectorMessenger;


class MyDetectorConstruction : public G4VUserDetectorConstruction {

public:

  MyDetectorConstruction();
  ~MyDetectorConstruction();

  G4VPhysicalVolume* Construct();

  void SetMagField( const G4double fieldValue );

private:
  
  G4VPhysicalVolume* fWorld;

  G4GDMLParser parser;

  G4FieldManager* fieldMgr;
  // Pointer to the field manager.

  G4UniformMagField* uniformMagField; 
  // Pointer to the uniform magnetic field.
  
  MyDetectorMessenger* detectorMessenger;
  // Pointer to the Messenger.

};

#endif

