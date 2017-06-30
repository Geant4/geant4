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
/// \file exoticphysics/monopole/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 104872 2017-06-23 14:19:16Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;
class G4MonopoleFieldSetup;
class G4GlobalMagFieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction();
  ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  // set geometry and field parameters
  void SetSizeX(G4double);
  void SetSizeYZ(G4double);              
  void SetMaterial(const G4String&);            
  void SetMagField   (G4double v) { fZMagFieldValue = v; }
  void SetMaxStepSize(G4double);
  void UpdateGeometry();

  virtual void ConstructSDandField();
          
  // access to geometry
  inline G4double     GetWorldSizeX()    {return fWorldSizeX;};
  inline G4double     GetAbsorSizeX()    {return fAbsorSizeX;};
  inline G4double     GetMaxStepSize()   {return fMaxStepSize;};
  inline const G4Material* GetAbsorMaterial() {return fAbsorMaterial;};

  G4MonopoleFieldSetup* GetMonopoleFieldSetup() const { return fMonFieldSetup; }
                           
private:

  void PrintParameters();
  
  G4double            fWorldSizeX;
  G4double            fWorldSizeYZ;
  G4Material*         fWorldMaterial;           
  G4double            fAbsorSizeX;
  G4double            fAbsorSizeYZ;
  G4double            fMaxStepSize;
  G4Material*         fAbsorMaterial;
  G4LogicalVolume*    fLogAbsor;

  G4MonopoleFieldSetup* fMonFieldSetup;
  G4double              fZMagFieldValue;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;
               
  DetectorMessenger*  fDetectorMessenger;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

