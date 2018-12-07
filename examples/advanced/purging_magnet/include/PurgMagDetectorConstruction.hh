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
// Code developed by:
//  S.Larsson
//
//    *****************************************
//    *                                       *
//    *    PurgMagDetectorConstruction.hh     *
//    *                                       *
//    *****************************************
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef PurgMagDetectorConstruction_h
#define PurgMagDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Cache.hh"
#include "G4MagneticField.hh"


class G4Box;
class G4Trd;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class PurgMagTabulatedField3D;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PurgMagDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  PurgMagDetectorConstruction();
  ~PurgMagDetectorConstruction();

public:  
     
  G4VPhysicalVolume* Construct();
  void ConstructSDandField();
     
public:

  void PrintDetectorParameters();
                    
  G4double GetWorldSizeXY()  const       {return WorldSizeXY;};
  G4double GetWorldSizeZ() const          {return WorldSizeZ;}; 

  G4double GetMeasureVolumeSizeXY() const {return MeasureVolumeSizeXY;}; 
  G4double GetMeasureVolumeSizeZ() const  {return MeasureVolumeSizeZ;};

  G4double GetGapSizeX1()            {return GapSizeX1;}; 
  G4double GetGapSizeX2()            {return GapSizeX2;}; 
  G4double GetGapSizeY1()            {return GapSizeY1;}; 
  G4double GetGapSizeY2()            {return GapSizeY2;}; 
  G4double GetGapSizeZ()             {return GapSizeZ;};

  G4Material* GetWorldMaterial()         {return WorldMaterial;};
  G4Material* GetGapMaterial()           {return GapMaterial;};
  
  const G4VPhysicalVolume* GetWorld() const          {return physiWorld;};           
  const G4VPhysicalVolume* GetMeasureVolume() const  {return physiMeasureVolume;};           
  const G4VPhysicalVolume* GetGap1() const           {return physiGap1;};
  const G4VPhysicalVolume* GetGap2() const      {return physiGap2;};

private:
   


  G4double           WorldSizeXY;
  G4double           WorldSizeZ;

  G4double           MeasureVolumeSizeXY;
  G4double           MeasureVolumeSizeZ;
  G4double           MeasureVolumePosition;

  G4double           GapSizeX1;
  G4double           GapSizeX2;
  G4double           GapSizeY1;
  G4double           GapSizeY2;
  G4double           GapSizeZ;
  G4double           Gap1PosX; 
  G4double           Gap1PosY; 
  G4double           Gap1PosZ; 
  G4double           Gap2PosX; 
  G4double           Gap2PosY; 
  G4double           Gap2PosZ; 

  G4double           SSD;
  G4double           zOffset;

  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;  
  G4Box*             solidWorld;
  

  G4VPhysicalVolume* physiGap1;
  G4LogicalVolume*   logicGap1;
  G4Trd*             solidGap1;


  G4VPhysicalVolume* physiGap2;
  G4LogicalVolume*   logicGap2;
  G4Trd*             solidGap2;

  G4VPhysicalVolume* physiMeasureVolume;
  G4LogicalVolume*   logicMeasureVolume;
  G4Box*             solidMeasureVolume;

  G4Material*        WorldMaterial;
  G4Material*        GapMaterial;

  G4Cache<G4MagneticField*> fField;  //pointer to the thread-local fields

private:
  


  void DefineMaterials();
  G4VPhysicalVolume* ConstructCalorimeter();     
};


#endif


