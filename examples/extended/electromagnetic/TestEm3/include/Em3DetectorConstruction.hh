// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3DetectorConstruction.hh,v 1.1 1999-10-11 16:55:48 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em3DetectorConstruction_h
#define Em3DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class Em3DetectorMessenger;
class Em3CalorimeterSD;

     const G4int MaxAbsor = 10;
     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em3DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Em3DetectorConstruction();
   ~Em3DetectorConstruction();

  public:
  
     void SetNbOfAbsor     (G4int);      
     void SetAbsorMaterial (G4int,G4String);     
     void SetAbsorThickness(G4int,G4double);     
     
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
     
     G4int       GetNbOfAbsor()             {return NbOfAbsor;}; 
     G4Material* GetAbsorMaterial(G4int i)  {return AbsorMaterial[i];};
     G4double    GetAbsorThickness(G4int i) {return AbsorThickness[i];};      
     
     const G4VPhysicalVolume* GetphysiWorld()        {return physiWorld;};           
     const G4VPhysicalVolume* GetAbsorber(G4int i)   {return physiAbsor[i];};
                 
  private:
     
     G4int              NbOfAbsor;     
     G4Material*        AbsorMaterial [MaxAbsor];
     G4double           AbsorThickness[MaxAbsor];
     
     G4int              NbOfLayers;
     G4double           LayerThickness;
          
     G4double           CalorSizeYZ;
     G4double           CalorThickness;
     
     G4Material*        defaultMaterial;
     G4double           WorldSizeYZ;
     G4double           WorldSizeX;
            
     G4Box*             solidWorld;       //pointer to the solid World 
     G4LogicalVolume*   logicWorld;       //pointer to the logical World
     G4VPhysicalVolume* physiWorld;       //pointer to the physical World

     G4Box*             solidCalor;       //pointer to the solid Calor 
     G4LogicalVolume*   logicCalor;       //pointer to the logical Calor
     G4VPhysicalVolume* physiCalor;       //pointer to the physical Calor
     
     G4Box*             solidLayer;       //pointer to the solid Layer 
     G4LogicalVolume*   logicLayer;       //pointer to the logical Layer
     G4VPhysicalVolume* physiLayer;       //pointer to the physical Layer
         
     G4Box*             solidAbsor[MaxAbsor]; //pointer to the solid Absorbers
     G4LogicalVolume*   logicAbsor[MaxAbsor]; //pointer to the logical Absorbers
     G4VPhysicalVolume* physiAbsor[MaxAbsor]; //pointer to the physical Absorbers
     
     G4UniformMagField* magField;         //pointer to the magnetic field
     
     Em3DetectorMessenger* detectorMessenger;  //pointer to the Messenger
     Em3CalorimeterSD* calorimeterSD;  //pointer to the sensitive detector
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Em3DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
     LayerThickness = 0.;
     for (G4int iAbs=0; iAbs<NbOfAbsor; iAbs++)
     LayerThickness += AbsorThickness[iAbs];
     CalorThickness = NbOfLayers*LayerThickness;
     
     WorldSizeX = 1.2*CalorThickness; WorldSizeYZ = 1.2*CalorSizeYZ;
}

#endif

