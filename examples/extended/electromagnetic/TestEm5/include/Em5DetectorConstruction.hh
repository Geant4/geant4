// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5DetectorConstruction.hh,v 1.1 1999-10-12 12:23:29 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5DetectorConstruction_h
#define Em5DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class Em5DetectorMessenger;
class Em5CalorimeterSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Em5DetectorConstruction();
   ~Em5DetectorConstruction();

  public:
     
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);     
     void SetAbsorberSizeYZ   (G4double);          
      
     void SetAbsorberXpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeX   (G4double);
     void SetWorldSizeYZ  (G4double);

     void SetMagField(G4double);
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintCalorParameters(); 
                    
     G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
     G4double    GetAbsorberThickness() {return AbsorberThickness;};      
     G4double    GetAbsorberSizeYZ()    {return AbsorberSizeYZ;};
     
     G4double    GetAbsorberXpos()      {return XposAbs;}; 
     G4double    GetxstartAbs()         {return xstartAbs;};
     G4double    GetxendAbs()           {return xendAbs;};
     
     G4Material* GetWorldMaterial()     {return WorldMaterial;};
     G4double    GetWorldSizeX()        {return WorldSizeX;}; 
     G4double    GetWorldSizeYZ()       {return WorldSizeYZ;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
                 
  private:
     
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;
     G4double           AbsorberSizeYZ;
 
     G4double           XposAbs;
     G4double           xstartAbs, xendAbs;
     
     G4Material*        WorldMaterial;
     G4double           WorldSizeX;     
     G4double           WorldSizeYZ;
     
     G4bool             defaultWorld;     
           
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     G4Box*             solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
     
     G4UniformMagField* magField;      //pointer to the magnetic field
     
     Em5DetectorMessenger* detectorMessenger;  //pointer to the Messenger
     Em5CalorimeterSD* calorimeterSD;  //pointer to the sensitive detector
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Em5DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
     
     xstartAbs = XposAbs-0.5*AbsorberThickness; 
     xendAbs   = XposAbs+0.5*AbsorberThickness;
     
     if (defaultWorld)
       {
        WorldSizeX = 1.5*AbsorberThickness; WorldSizeYZ= 1.2*AbsorberSizeYZ;
       } 	
}

#endif

