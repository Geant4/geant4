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
// $Id: DetectorConstruction.hh,v 1.1 2005/06/03 15:19:49 maire Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
       
     void SetSizeX      (G4double);
     void SetSizeYZ     (G4double);           
     void SetMaterial   (G4String);
     void SetNbOfLayers (G4int nb);                    
     void SetMagField   (G4double);
          
     G4VPhysicalVolume* Construct();
     void               UpdateGeometry();
     
  public:  
                    
     G4double     GetAbsorSizeX()    {return absorSizeX;};
     G4double     GetAbsorSizeYZ()   {return absorSizeYZ;};           
     G4Material*  GetAbsorMaterial() {return absorMaterial;};
     G4int        GetNbOfLayers()    {return nbOfLayers;};   
     
     void         PrintParameters();
                       
  private:

     G4double            absorSizeX;
     G4double            absorSizeYZ;     
     G4Material*         absorMaterial;
     G4int               nbOfLayers;
     G4double            layerThickness;
     G4UniformMagField*  magField;
     G4VPhysicalVolume*  pAbsor;

     DetectorMessenger* detectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

