//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst14DetectorConstruction.hh,v 1.6 2003-01-31 08:06:11 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst14DetectorConstruction_h
#define Tst14DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class Tst14DetectorMessenger;
class Tst14CalorimeterSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst14DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Tst14DetectorConstruction();
   ~Tst14DetectorConstruction();

  public:
     
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);     
     void SetAbsorberRadius(G4double);          
      
     void SetAbsorberZpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);

     void SetMagField(G4double);
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintCalorParameters(); 
                    
     G4Material* GetWorldMaterial()    {return WorldMaterial;};
     G4double GetWorldSizeZ()          {return WorldSizeZ;}; 
     G4double GetWorldSizeR()          {return WorldSizeR;};
     
     G4double GetAbsorberZpos()        {return zAbsorber;}; 
     G4double GetzstartAbs()           {return zstartAbs;};
     G4double GetzendAbs()             {return zendAbs;};

     G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
     G4double    GetAbsorberThickness() {return AbsorberThickness;};      
     G4double GetAbsorberRadius()       {return AbsorberRadius;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
                 
  private:
     
     G4bool             worldchanged;
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;
     G4double           AbsorberRadius;
 
     G4double           zAbsorber ;
     G4double           zstartAbs , zendAbs ;
     
     G4Material*        WorldMaterial;
     G4double           WorldSizeR;
     G4double           WorldSizeZ;
            
     G4Tubs*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     G4Tubs*             solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
     
     G4UniformMagField* magField;      //pointer to the magnetic field
     
     Tst14DetectorMessenger* detectorMessenger;  //pointer to the Messenger
     Tst14CalorimeterSD* calorimeterSD;  //pointer to the sensitive detector
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};



#endif

