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
<<<<<<< HEAD
// $Id: DetectorConstruction.hh 68036 2013-03-13 14:13:45Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

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

  virtual G4VPhysicalVolume* Construct();

  // set geometry and field parameters
  void SetSizeX(G4double);
  void SetSizeYZ(G4double);              
  void SetMaterial(const G4String&);            
<<<<<<< HEAD
  void SetMagField(G4double);
  void SetMaxStepSize(G4double);
  void UpdateGeometry();
=======
  void SetMaxStepSize(G4double);

  virtual void ConstructSDandField();
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
          
  // access to geometry
  inline G4double     GetWorldSizeX()    {return fWorldSizeX;};
  inline G4double     GetAbsorSizeX()    {return fAbsorSizeX;};
  inline G4double     GetMaxStepSize()   {return fMaxStepSize;};
  inline const G4Material* GetAbsorMaterial() {return fAbsorMaterial;};
<<<<<<< HEAD
=======

  G4MonopoleFieldSetup* GetMonopoleFieldSetup() const { return fMonFieldSetup.Get(); }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
                           
private:

  void PrintParameters();
  
  G4double            fWorldSizeX;
  G4double            fWorldSizeYZ;
  G4Material*         fWorldMaterial;           
  G4double            fAbsorSizeX;
  G4double            fAbsorSizeYZ;
  G4double            fMaxStepSize;
  G4Material*         fAbsorMaterial;

<<<<<<< HEAD
  G4UniformMagField*    fMagField;
  G4MonopoleFieldSetup* fMonFieldSetup;

  G4LogicalVolume*    fLogAbsor;
=======
  G4Cache<G4MonopoleFieldSetup*> fMonFieldSetup;
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
               
  DetectorMessenger*  fDetectorMessenger;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

