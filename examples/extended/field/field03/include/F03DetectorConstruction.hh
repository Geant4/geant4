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
// $Id: F03DetectorConstruction.hh,v 1.6 2006-06-29 17:18:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef F03DetectorConstruction_h
#define F03DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class F03DetectorMessenger;
class F03CalorimeterSD;
class F03FieldSetup;


class F03DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    F03DetectorConstruction();
   ~F03DetectorConstruction();

  public:
     
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);     
     void SetAbsorberRadius(G4double);          
      
     void SetAbsorberZpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);
     
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
     
     G4Tubs*            solidWorld;    // pointer to the solid World 
     G4LogicalVolume*   logicWorld;    // pointer to the logical World
     G4VPhysicalVolume* physiWorld;    // pointer to the physical World

     G4Tubs*            solidAbsorber; // pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber; // pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; // pointer to the physical Absorber
     
     G4UniformMagField* magField;      // pointer to the magnetic field
     F03FieldSetup*        fEmFieldSetup;     
     F03DetectorMessenger* detectorMessenger;  // pointer to the Messenger
     F03CalorimeterSD* calorimeterSD;  // pointer to the sensitive detector
      
  // TR radiator volumes and dimensions
          
     G4Tubs*            fSolidRadSlice;   // pointer to the solid  z-slice 
     G4LogicalVolume*   fLogicRadSlice;   // pointer to the logical z-slide
     G4VPhysicalVolume* fPhysicRadSlice;  // pointer to the physical z-slide

     G4Tubs*            solidRadiator;
     G4LogicalVolume*   logicRadiator; 
     G4VPhysicalVolume* physiRadiator;

     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;
     G4double           AbsorberRadius;
     G4Material*        fRadiatorMat;     // pointer to the TR radiator material
     G4bool             worldchanged;

     G4double           zAbsorber ;
     G4double           zstartAbs , zendAbs ;
     
     G4Material*        WorldMaterial;
     G4double           WorldSizeR;
     G4double           WorldSizeZ;

     G4double fRadThickness ;
     G4double fGasGap       ;

     G4int fFoilNumber ;

     G4double fDetGap       ;

     G4double fStartR       ;
     G4double fStartZ       ;

  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void F03DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
     
     zstartAbs = zAbsorber-0.5*AbsorberThickness; 
     zendAbs   = zAbsorber+0.5*AbsorberThickness; 

}

#endif
