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
// $Id: F02DetectorConstruction.hh,v 1.2 2001-07-11 09:58:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef F02DetectorConstruction_h
#define F02DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class F02DetectorMessenger;
class F02CalorimeterSD;
class F02ElectroMagneticField;


class F02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    F02DetectorConstruction();
   ~F02DetectorConstruction();

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
     G4LogicalVolume* GetLogicalAbsorber()    {return logicAbsorber;};
                 
  private:
     
     G4bool             worldchanged;
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;
     G4double           AbsorberRadius;

  G4Material*        fWindowMat ;
  G4double           fWindowThick ;

  G4Material*        fElectrodeMat ;
  G4double           fElectrodeThick ;

  G4Material*        fGapMat ;
  G4double           fGapThick ;

 
     G4double           zAbsorber ;
     G4double           zstartAbs , zendAbs ;
     
     G4Material*        WorldMaterial;
     G4double           WorldSizeR;
     G4double           WorldSizeZ;
            
     G4Tubs*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

  // TR radiator volumes and dimensions
          
     G4Tubs*            fSolidRadSlice;   // pointer to the solid  z-slice 
     G4LogicalVolume*   fLogicRadSlice;   // pointer to the logical z-slide
     G4VPhysicalVolume* fPhysicRadSlice;  // pointer to the physical z-slide

     G4Tubs*            fSolidRadRing;    // pointer to the solid  R-slice 
     G4LogicalVolume*   fLogicRadRing;    // pointer to the logical R-slide
     G4VPhysicalVolume* fPhysicRadRing;   // pointer to the physical R-slide

     G4Material* fRadiatorMat;        //pointer to the TR radiator material

     G4double fRadThickness ;
     G4double fGasGap       ;

     G4int fFoilNumber ;

     G4double fDetThickness ;
     G4double fDetLength    ;
     G4double fDetGap       ;

     G4double fStartR       ;
     G4double fStartZ       ;

     G4int fModuleNumber ;   // the number of Rad-Det modules

     G4Tubs*             solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
     
     G4UniformMagField* magField;      //pointer to the magnetic field
     F02ElectroMagneticField* fEmField;     
     F02DetectorMessenger* detectorMessenger;  //pointer to the Messenger
     F02CalorimeterSD* calorimeterSD;  //pointer to the sensitive detector
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void F02DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
     if(!worldchanged)
     {
       //  WorldSizeR=2.*AbsorberRadius ;
       //  WorldSizeZ=2.*AbsorberThickness ;
     }
     
     zstartAbs = zAbsorber-0.5*AbsorberThickness; 
     zendAbs   = zAbsorber+0.5*AbsorberThickness; 

}

#endif

