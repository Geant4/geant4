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
/// \file field/field02/include/F02DetectorConstruction.hh
/// \brief Definition of the F02DetectorConstruction class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F02DetectorConstruction_h
#define F02DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"

#include "CLHEP/Units/SystemOfUnits.h"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;

class G4Material;
class G4UniformMagField;

class F02DetectorMessenger;
class F02CalorimeterSD;
class F02ElectricFieldSetup;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    F02DetectorConstruction();
    ~F02DetectorConstruction() override;

  public:

     void SetAbsorberMaterial (G4String);
     void SetAbsorberThickness(G4double);
     void SetAbsorberRadius(G4double);

     void SetAbsorberZpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);

     G4VPhysicalVolume* Construct() override;
     void ConstructSDandField() override;

  public:

     void PrintCalorParameters();

     G4Material* GetWorldMaterial()    {return fWorldMaterial;}
     G4double GetWorldSizeZ()          {return fWorldSizeZ;}
     G4double GetWorldSizeR()          {return fWorldSizeR;}

     G4double GetAbsorberZpos()        {return fZAbsorber;}
     G4double GetZStartAbs()           {return fZStartAbs;}
     G4double GetZEndAbs()             {return fZEndAbs;}

     G4Material* GetAbsorberMaterial() {return fAbsorberMaterial;}
     G4double    GetAbsorberThickness(){return fAbsorberThickness;}
     G4double    GetAbsorberRadius()   {return fAbsorberRadius;}

     const G4VPhysicalVolume* GetPhysiWorld() {return fPhysiWorld;}
     const G4VPhysicalVolume* GetAbsorber()   {return fPhysiAbsorber;}
     G4LogicalVolume* GetLogicalAbsorber()    {return fLogicAbsorber;}

  private:

     F02DetectorMessenger* fDetectorMessenger = nullptr;  // pointer -> Messenger
     G4Cache<F02CalorimeterSD*> fCalorimeterSD = nullptr; // pointer -> sensitive detector
     G4Cache<F02ElectricFieldSetup*> fEmFieldSetup = nullptr;

     G4Tubs*            fSolidWorld = nullptr;     // pointer to the solid World
     G4LogicalVolume*   fLogicWorld = nullptr;     // pointer to the logical World
     G4VPhysicalVolume* fPhysiWorld = nullptr;     // pointer to the physical World

     G4Tubs*            fSolidAbsorber = nullptr;  // pointer to the solid Absorber
     G4LogicalVolume*   fLogicAbsorber = nullptr;  // pointer to the logical Absorber
     G4VPhysicalVolume* fPhysiAbsorber = nullptr;  // pointer to the physical Absorber

     G4Material*        fAbsorberMaterial = nullptr;
     G4double           fAbsorberThickness  = 4. * CLHEP::cm;
     G4double           fAbsorberRadius  = 10. * CLHEP::cm;
     G4bool             fWorldChanged;

     G4double           fZAbsorber = 36. * CLHEP::cm;
     G4double           fZStartAbs = 0.;
     G4double           fZEndAbs = 0.;

     G4Material*        fWorldMaterial = nullptr;
     G4double           fWorldSizeR = 20. * CLHEP::cm;
     G4double           fWorldSizeZ = 80. * CLHEP::cm;

  private:

     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void F02DetectorConstruction::ComputeCalorParameters()
{
     // Compute derived parameters of the calorimeter

     fZStartAbs = fZAbsorber-0.5*fAbsorberThickness;
     fZEndAbs   = fZAbsorber+0.5*fAbsorberThickness;
}

#endif
