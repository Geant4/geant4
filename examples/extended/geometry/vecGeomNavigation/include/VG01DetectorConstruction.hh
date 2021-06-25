
#ifndef VG01DetectorConstruction_h
#define VG01DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "G4GDMLParser.hh"
#include "G4String.hh"

class G4VPhysicalVolume;
class G4FieldManager;
class G4UniformMagField;
class VG01DetectorMessenger;

class VG01DetectorConstruction : public G4VUserDetectorConstruction {

public:

  VG01DetectorConstruction();
  ~VG01DetectorConstruction();

  G4VPhysicalVolume* Construct() override;

  void SetGDMLFileName ( const G4String& gdmlfile ) { fGDMLFileName = gdmlfile; }
  void SetUseVecGeom (bool b) { fUseVecGeom = b; }
  void SetMagFieldValue(const G4double fieldValue ) { fglobFieldValue = fieldValue; }

  static G4double GetFieldValue() { return fglobFieldValue; }

private:
  void CreateMagFieldAndIntegrator();

private:
  // this static member is for the print out
  static G4double        fglobFieldValue;

  G4String               fGDMLFileName;
  G4GDMLParser           fParser;
  G4VPhysicalVolume*     fWorld;
  G4FieldManager*        fFieldMgr;
  G4UniformMagField*     fUniformMagField;
  bool                   fUseVecGeom = true;
  VG01DetectorMessenger* fDetectorMessenger= nullptr;
};

#endif
