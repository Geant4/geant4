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
// $Id: DetectorConstruction.hh,v 1.1 2003/04/22 16:25:04 maire Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

      const G4int MaxTally = 20;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
       
     void SetSizeX    (G4double);
     void SetSizeYZ   (G4double);              
     void SetMaterial (G4String);            
     void SetMagField (G4double);
     
     void SetTallySize     (G4ThreeVector);
     void SetTallyMaterial (G4String);
     void SetTallyPosition (G4ThreeVector);
     
     G4VPhysicalVolume* Construct();
     void               UpdateGeometry();
     
  public:  
                    
     G4double     GetWorldSizeX()    {return worldSizeX;};
     G4double     GetWorldSizeYZ()   {return worldSizeYZ;};
     G4Material*  GetWorldMaterial() {return worldMaterial;};     
     G4double     GetAbsorSizeX()    {return absorSizeX;};
     G4double     GetAbsorSizeYZ()   {return absorSizeYZ;};           
     G4Material*  GetAbsorMaterial() {return absorMaterial;};
     
     G4LogicalVolume* GetLogicalTally() {return lTally;}
     G4double         GetTallyMass()    {return tallyMass;};
     G4int            GetTallyNumber()  {return tallyNumber;};
     
     void         PrintParameters();
                       
  private:
  
     G4double            worldSizeX;
     G4double            worldSizeYZ;
     G4Material*         worldMaterial;           
     G4double            absorSizeX;
     G4double            absorSizeYZ;     
     G4Material*         absorMaterial;
     G4UniformMagField*  magField;
     G4LogicalVolume*    lAbsor;
               
     G4ThreeVector       tallySize;
     G4Material*         tallyMaterial;
     G4double            tallyMass;
     G4int               tallyNumber;     
     G4ThreeVector*      tallyPosition;
     G4LogicalVolume*    lTally;
     
     DetectorMessenger* detectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

