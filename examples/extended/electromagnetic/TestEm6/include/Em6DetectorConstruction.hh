//  Em6DetectorConstruction.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em6DetectorConstruction_h
#define Em6DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"

class G4Box;
class G4LogicalVolume;
class G4UniformMagField;
class Em6DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    Em6DetectorConstruction();
   ~Em6DetectorConstruction();

  public:

     void SetMaterial(G4String);
     void SetLBining (G4ThreeVector);
     void SetMagField(G4double);

     G4VPhysicalVolume* Construct();

     void UpdateGeometry();

  public:

     const
     G4VPhysicalVolume* GetAbsorber() {return physiAbsorber;};
     G4Material*    GetMaterial() {return defaultMaterial;};

     G4int    GetnLtot()      {return nLtot;};
     G4double GetdLradl()     {return dLradl;};
     G4double GetfullLength() {return AbsorberLength;};

  private:

     G4int    nLtot;          // nb of bins: longitudinal
     G4double dLradl;         // bin thickness (in radl unit)

     G4Material* defaultMaterial;     //pointer to the material
     G4UniformMagField* magField;     //pointer to the mag field

     G4double AbsorberLength;       //full length of the Absorber
     G4double AbsorberSizeXY;       //full x,y size of the Absorber

     G4Box*             solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber

     G4Box*            solidLayer;   //pointer to the solid  L-Layer
     G4LogicalVolume*   logicLayer;  //pointer to the logical L-slide
     G4VPhysicalVolume* physiLayer;  //pointer to the physical L-slide

     Em6DetectorMessenger* detectorMessenger;  //pointer to the Messenger

  private:

     void DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

