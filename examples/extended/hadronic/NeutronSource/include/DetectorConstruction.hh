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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
    virtual G4VPhysicalVolume* Construct();

    G4Material* 
    MaterialWithSingleIsotope(G4String, G4String, G4double, G4int, G4int);
         
    void SetAbsorRadius   (G4double);
    void SetAbsorLength   (G4double);
    void SetAbsorMaterial (G4String);

    void SetContainThickness (G4double);
    void SetContainMaterial (G4String);

  public:

     G4double           GetAbsorRadius()     {return fAbsorRadius;};
     G4double           GetAbsorLength()     {return fAbsorLength;};
     G4Material*        GetAbsorMaterial()   {return fAbsorMaterial;};

     G4double           GetContainThickness()  {return fContainThickness;};
     G4Material*        GetContainMaterial()   {return fContainMaterial;};

     void               PrintParameters();

  private:

     G4double           fAbsorRadius, fAbsorLength;
     G4Material*        fAbsorMaterial;
     G4LogicalVolume*   fLAbsor;
     
     G4double           fContainThickness;
     G4Material*        fContainMaterial;
     G4LogicalVolume*   fLContain;

     G4Material*        fWorldMaterial;
     G4VPhysicalVolume* fPWorld;
     
     DetectorMessenger* fDetectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

