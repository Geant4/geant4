// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03DetectorConstruction.hh,v 1.1 1999-01-07 16:05:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ExN03DetectorConstruction_h
#define ExN03DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class ExN03DetectorMessenger;
class ExN03CalorimeterSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ExN03DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    ExN03DetectorConstruction();
   ~ExN03DetectorConstruction();

  public:
     
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);     

     void SetGapMaterial (G4String);     
     void SetGapThickness(G4double);
     
     void SetCalorSizeYZ(G4double);          
     void SetNbOfLayers (G4int);   
      
     void SetMagField(G4double);
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintCalorParameters(); 
                    
     G4double GetWorldSizeX()           {return WorldSizeX;}; 
     G4double GetWorldSizeYZ()          {return WorldSizeYZ;};
     
     G4double GetCalorThickness()       {return CalorThickness;}; 
     G4double GetCalorSizeYZ()          {return CalorSizeYZ;};
      
     G4int GetNbOfLayers()              {return NbOfLayers;}; 
     
     G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
     G4double    GetAbsorberThickness() {return AbsorberThickness;};      
     
     G4Material* GetGapMaterial()       {return GapMaterial;};
     G4double    GetGapThickness()      {return GapThickness;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
     const G4VPhysicalVolume* GetGap()        {return physiGap;};
                 
  private:
     
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;
     
     G4Material*        GapMaterial;
     G4double           GapThickness;
     
     G4int              NbOfLayers;
     G4double           LayerThickness;
          
     G4double           CalorSizeYZ;
     G4double           CalorThickness;
     
     G4Material*        defaultMaterial;
     G4double           WorldSizeYZ;
     G4double           WorldSizeX;
            
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     G4Box*             solidCalor;    //pointer to the solid Calor 
     G4LogicalVolume*   logicCalor;    //pointer to the logical Calor
     G4VPhysicalVolume* physiCalor;    //pointer to the physical Calor
     
     G4Box*             solidLayer;    //pointer to the solid Layer 
     G4LogicalVolume*   logicLayer;    //pointer to the logical Layer
     G4VPhysicalVolume* physiLayer;    //pointer to the physical Layer
         
     G4Box*             solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
     
     G4Box*             solidGap;      //pointer to the solid Gap
     G4LogicalVolume*   logicGap;      //pointer to the logical Gap
     G4VPhysicalVolume* physiGap;      //pointer to the physical Gap
     
     G4UniformMagField* magField;      //pointer to the magnetic field
     
     ExN03DetectorMessenger* detectorMessenger;  //pointer to the Messenger
     ExN03CalorimeterSD* calorimeterSD;  //pointer to the sensitive detector
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void ExN03DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
     LayerThickness = AbsorberThickness + GapThickness;
     CalorThickness = NbOfLayers*LayerThickness;
     
     WorldSizeX = 1.2*CalorThickness; WorldSizeYZ = 1.2*CalorSizeYZ;
}

#endif

