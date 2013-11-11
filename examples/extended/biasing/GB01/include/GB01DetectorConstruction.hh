#ifndef GB01DetectorConstruction_h
#define GB01DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GB01DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  GB01DetectorConstruction();
  ~GB01DetectorConstruction();
  
public:
  
  virtual G4VPhysicalVolume* Construct();

};

#endif

