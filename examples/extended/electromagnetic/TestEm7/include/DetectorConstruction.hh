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
/// \file electromagnetic/TestEm7/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 66241 2012-12-13 18:34:42Z gunter $
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

      const G4int MaxTally = 21;                // 0 + 20

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
     
     void SetTallyNumber   (G4int);     
     void SetTallySize     (G4int, G4ThreeVector);
     void SetTallyMaterial (G4int, G4String);
     void SetTallyPosition (G4int, G4ThreeVector);

     virtual     
     G4VPhysicalVolume* Construct();
     void               UpdateGeometry();
     
  public:  
                    
     G4double     GetWorldSizeX()    {return fWorldSizeX;};
     G4double     GetWorldSizeYZ()   {return fWorldSizeYZ;};
     G4Material*  GetWorldMaterial() {return fWorldMaterial;};     
     G4double     GetAbsorSizeX()    {return fAbsorSizeX;};
     G4double     GetAbsorSizeYZ()   {return fAbsorSizeYZ;};           
     G4Material*  GetAbsorMaterial() {return fAbsorMaterial;};
     
     G4int            GetTallyNumber()         {return fTallyNumber;};
     G4double         GetTallyMass(G4int n)    {return fTallyMass[n];};          
     G4LogicalVolume* GetLogicalTally(G4int n) {return fLTally[n];}
     
     void         PrintParameters();
                       
  private:
  
     G4double            fWorldSizeX;
     G4double            fWorldSizeYZ;
     G4Material*         fWorldMaterial;           
     G4double            fAbsorSizeX;
     G4double            fAbsorSizeYZ;     
     G4Material*         fAbsorMaterial;
     G4UniformMagField*  fMagField;
     G4LogicalVolume*    fLAbsor;
     
     G4int               fTallyNumber;                   
     G4ThreeVector       fTallySize[MaxTally];
     G4Material*         fTallyMaterial[MaxTally];
     G4double            fTallyMass[MaxTally]; 
     G4ThreeVector       fTallyPosition[MaxTally];
     G4LogicalVolume*    fLTally[MaxTally];
     
     DetectorMessenger*  fDetectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

