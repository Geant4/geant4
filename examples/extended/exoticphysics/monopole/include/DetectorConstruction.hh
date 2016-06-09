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
// $Id: DetectorConstruction.hh,v 1.3 2010-06-04 19:03:36 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
class G4MonopoleFieldSetup;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction();
  ~DetectorConstruction();

  void SetSizeX    (G4double);
  void SetSizeYZ   (G4double);              
  void SetMaterial (G4String);            
  void SetMagField (G4double);
     
  void SetMaxStepSize   (G4double);
     
  G4VPhysicalVolume* Construct();
  void               UpdateGeometry();
     
  G4double     GetWorldSizeX()    {return worldSizeX;};
  G4double     GetWorldSizeYZ()   {return worldSizeYZ;};
  G4Material*  GetWorldMaterial() {return worldMaterial;};     
  G4double     GetAbsorSizeX()    {return absorSizeX;};
  G4double     GetAbsorSizeYZ()   {return absorSizeYZ;};
  G4double     GetMaxStepSize()   {return maxStepSize;};
  G4Material*  GetAbsorMaterial() {return absorMaterial;};
  
  void         PrintParameters();
                       
private:

  void               DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();     
  
  G4double            worldSizeX;
  G4double            worldSizeYZ;
  G4Material*         worldMaterial;           
  G4double            absorSizeX;
  G4double            absorSizeYZ;
  G4double	      maxStepSize;
  G4Material*         absorMaterial;

  G4UniformMagField*  magField;
  G4MonopoleFieldSetup* fMFieldSetup;

  G4LogicalVolume*    lAbsor;
               
  DetectorMessenger* detectorMessenger;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

