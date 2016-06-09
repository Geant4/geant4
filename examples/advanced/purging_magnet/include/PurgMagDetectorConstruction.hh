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
// Code developed by:
//  S.Larsson
//
//    *****************************************
//    *                                       *
//    *    PurgMagDetectorConstruction.hh     *
//    *                                       *
//    *****************************************
//
// $Id: PurgMagDetectorConstruction.hh,v 1.2 2004/06/18 09:17:44 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef PurgMagDetectorConstruction_h
#define PurgMagDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"


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
     
public:

  void PrintDetectorParameters();
                    
  G4double GetWorldSizeXY()          {return WorldSizeXY;};
  G4double GetWorldSizeZ()           {return WorldSizeZ;}; 

  G4double GetMeasureVolumeSizeXY()  {return MeasureVolumeSizeXY;}; 
  G4double GetMeasureVolumeSizeZ()   {return MeasureVolumeSizeZ;};

  G4double GetGapSizeX1()            {return GapSizeX1;}; 
  G4double GetGapSizeX2()            {return GapSizeX2;}; 
  G4double GetGapSizeY1()            {return GapSizeY1;}; 
  G4double GetGapSizeY2()            {return GapSizeY2;}; 
  G4double GetGapSizeZ()             {return GapSizeZ;};

  G4Material* GetWorldMaterial()         {return WorldMaterial;};
  G4Material* GetGapMaterial()           {return GapMaterial;};
  
  const G4VPhysicalVolume* GetWorld()           {return physiWorld;};           
  const G4VPhysicalVolume* GetMeasureVolume()   {return physiMeasureVolume;};           
  const G4VPhysicalVolume* GetGap1()            {return physiGap1;};
  const G4VPhysicalVolume* GetGap2()            {return physiGap2;};

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



private:
    
  void DefineMaterials();
  G4VPhysicalVolume* ConstructCalorimeter();     
};


#endif


