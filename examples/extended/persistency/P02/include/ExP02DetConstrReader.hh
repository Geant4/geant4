//
#ifndef ExP02DetConstrReader_h
#define ExP02DetConstrReader_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExP02DetConstrReader : public G4VUserDetectorConstruction
{
  public:
  
     ExP02DetConstrReader();
    ~ExP02DetConstrReader();

  public:
  
     G4VPhysicalVolume* Construct();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
