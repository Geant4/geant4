#ifndef hTestDetectorConstruction_h
#define hTestDetectorConstruction_h 1

// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestDetectorConstruction -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class hTestDetectorMessenger;
class hTestCalorimeterSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    hTestDetectorConstruction();
   ~hTestDetectorConstruction();

  public:
     
    void SetAbsorberMaterial (G4String);     
    void SetNumberOfAbsorbers(G4int);     
    void SetAbsorberThickness(G4double);     
    void SetAbsorberSizeXY   (G4double);          
      
    void SetWorldMaterial(G4String);
    void SetWorldSizeZ   (G4double);

    void SetMagField(G4double,G4int);
    void SetVerbose(G4int val) {verbose = val;};
    void SetHistoName(G4String name) {histoName = name;};
     
    G4VPhysicalVolume* Construct();

    void UpdateGeometry();
     
    void PrintGeomParameters(); 
                    
    G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
    G4double    GetAbsorberThickness() {return AbsorberThickness;};      
    G4double    GetAbsorberSizeXY()    {return SizeXY;};
          
    G4Material* GetWorldMaterial()     {return WorldMaterial;};
    G4double    GetWorldSizeZ()        {return WorldSizeZ;}; 
    G4double    GetWorldSizeXY()       {return SizeXY;};
     
    const G4VPhysicalVolume* GetPhysWorld() {return physWorld;};           
    const G4LogicalVolume*   GetAbsorber()  {return logicAbs;};

    G4String GetHistoName() const {return histoName;};
                 
  private:

    // Methods

     void DefineMaterials();
     void ComputeGeomParameters();
     G4VPhysicalVolume* ConstructGeometry();     
     void GeometryIsChanged();     
     void MaterialIsChanged();     

    // Members
     
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;  // 
     G4double           SizeXY;

     G4int              NumberOfAbsorbers;

     G4double           ZposAbs;
     
     G4Material*        WorldMaterial;
     G4double           WorldSizeZ;     
                
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physWorld;     //pointer to the physical World

     G4Box*             solidAbs;      //pointer to the solid Absorber
     G4LogicalVolume*   logicAbs;      //pointer to the logical Absorber
     G4VPhysicalVolume* physAbs;       //pointer to the physical Absorber
     
     G4UniformMagField* magField;      //pointer to the magnetic field
     
     hTestDetectorMessenger* detectorMessenger;  //pointer to the Messenger
     hTestCalorimeterSD* calorimeterSD;  //pointer to the sensitive detector

     G4int verbose;      
     G4String histoName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

