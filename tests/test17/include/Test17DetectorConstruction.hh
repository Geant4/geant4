//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------
//      GEANT 4 class example
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- Test17DetectorConstruction -------
//                by Vladimir Ivanchenko, 23 July 1999 
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17DetectorConstruction_h
#define Test17DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class Test17DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Test17DetectorConstruction();
   ~Test17DetectorConstruction();

  public:
     
     void SetAbsorberMaterial (const G4String&);     
     void SetNumberOfAbsorbers (G4int);     
     void SetAbsorberThickness(G4double);     
     void SetAbsorberSizeYZ   (G4double);          
      
     void SetAbsorberXpos(G4double);

     void SetWorldMaterial(const G4String&);
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
     const G4LogicalVolume* GetAbsorber()   {return logicAbsorber;};
                 
  private:
     
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;
     G4double           AbsorberSizeYZ;

     G4int              NumberOfAbsorbers;

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
     
     Test17DetectorMessenger* detectorMessenger;  //pointer to the Messenger
  //Test17CalorimeterSD* calorimeterSD;  //pointer to the sensitive detector
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Test17DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the 1st absorber 
     
     xstartAbs = XposAbs-0.5*AbsorberThickness; 
     xendAbs   = XposAbs+0.5*AbsorberThickness;
     
     if (defaultWorld)
       {
        WorldSizeX  = 2.1*AbsorberThickness*NumberOfAbsorbers; 
        WorldSizeYZ = 1.2*AbsorberSizeYZ;
       } 	
}

#endif

 


