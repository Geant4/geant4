#ifndef GB02DetectorConstruction_h
#define GB02DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GB02DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  GB02DetectorConstruction();
  ~GB02DetectorConstruction();
  
public:
  
  virtual G4VPhysicalVolume* Construct();
  virtual void     ConstructSDandField();
  
};

#endif

