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
/// \file field/field03/include/F03DetectorConstruction.hh
/// \brief Definition of the F03DetectorConstruction class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F03DetectorConstruction_h
#define F03DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"

#include "CLHEP/Units/SystemOfUnits.h"

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;

class G4Material;
class G4UniformMagField;

class F03DetectorMessenger;
class F03CalorimeterSD;
class F03FieldSetup;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F03DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    F03DetectorConstruction();
    ~F03DetectorConstruction() override;

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
     F03DetectorMessenger* fDetectorMessenger = nullptr;  // pointer to the Messenger
     G4Cache<F03CalorimeterSD*> fCalorimeterSD = nullptr; // pointer to the sensitive det.
     G4Cache<F03FieldSetup*>    fEmFieldSetup = nullptr;

     G4Tubs*            fSolidWorld = nullptr;     // pointer to the solid World
     G4LogicalVolume*   fLogicWorld = nullptr;     // pointer to the logical World
     G4VPhysicalVolume* fPhysiWorld = nullptr;     // pointer to the physical World

     G4Tubs*            fSolidAbsorber = nullptr;  // pointer to the solid Absorber
     G4LogicalVolume*   fLogicAbsorber = nullptr;  // pointer to the logical Absorber
     G4VPhysicalVolume* fPhysiAbsorber = nullptr;  // pointer to the physical Absorber

     G4Tubs*            fSolidRadSlice = nullptr;  // pointer to the solid  z-slice
     G4LogicalVolume*   fLogicRadSlice = nullptr;  // pointer to the logical z-slide
     G4VPhysicalVolume* fPhysiRadSlice = nullptr;  // pointer to the physical z-slide

     G4Tubs*            fSolidRadiator = nullptr;
     G4LogicalVolume*   fLogicRadiator = nullptr;
     G4VPhysicalVolume* fPhysiRadiator = nullptr;

     G4Material*        fWorldMaterial = nullptr;
     G4Material*        fAbsorberMaterial = nullptr;
     G4Material*        fRadiatorMat = nullptr;    // pointer to the TR radiator material

     G4double           fWorldSizeR = 22000. * CLHEP::mm;
     G4double           fWorldSizeZ = 44000. * CLHEP::mm;

     G4double           fAbsorberThickness =  1. * CLHEP::mm;
     G4double           fAbsorberRadius = 20000. * CLHEP::mm;

     G4double           fZAbsorber = 21990. * CLHEP::mm;
     G4double           fZStartAbs = 0.;
     G4double           fZEndAbs =   0.;

     G4double           fRadThickness = 100. * CLHEP::mm;
     G4double           fGasGap =       100. * CLHEP::mm;
     G4double           fDetGap =         1. * CLHEP::mm;

     G4int              fFoilNumber = 2;

  private:

     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void F03DetectorConstruction::ComputeCalorParameters()
{
     // Compute derived parameters of the calorimeter

     fZStartAbs = fZAbsorber-0.5*fAbsorberThickness;
     fZEndAbs   = fZAbsorber+0.5*fAbsorberThickness;
}

#endif
