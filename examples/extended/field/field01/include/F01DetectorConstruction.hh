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
/// \file field/field01/include/F01DetectorConstruction.hh
/// \brief Definition of the F01DetectorConstruction class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F01DetectorConstruction_h
#define F01DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;

class G4Material;
class G4UniformMagField;

class F01DetectorMessenger;
class F01CalorimeterSD;
class F01FieldSetup;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    F01DetectorConstruction();
    virtual ~F01DetectorConstruction();

  public:

     void SetAbsorberMaterial (G4String);
     void SetAbsorberThickness(G4double);
     void SetAbsorberRadius(G4double);

     void SetAbsorberZpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);

     virtual G4VPhysicalVolume* Construct();
     virtual void ConstructSDandField();

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

     void    SetUseFSALstepper( G4bool val ) { fUseFSALstepper = val; }
     G4bool  AreUsingFSALstepper() { return fUseFSALstepper; }
   
  private:

     F01DetectorMessenger* fDetectorMessenger;  // pointer -> Messenger
     G4Cache<F01CalorimeterSD*> fCalorimeterSD; // pointer -> sensitive detector
     G4Cache<F01FieldSetup*>    fEmFieldSetup;

     G4Tubs*            fSolidWorld;     // pointer to the solid World
     G4LogicalVolume*   fLogicWorld;     // pointer to the logical World
     G4VPhysicalVolume* fPhysiWorld;     // pointer to the physical World

     G4Tubs*            fSolidAbsorber;  // pointer to the solid Absorber
     G4LogicalVolume*   fLogicAbsorber;  // pointer to the logical Absorber
     G4VPhysicalVolume* fPhysiAbsorber;  // pointer to the physical Absorber
 
     G4Material*        fAbsorberMaterial;
     G4double           fAbsorberThickness;
     G4double           fAbsorberRadius;

     G4double           fZAbsorber;
     G4double           fZStartAbs, fZEndAbs;

     G4Material*        fWorldMaterial;
     G4double           fWorldSizeR;
     G4double           fWorldSizeZ;

     G4bool             fUseFSALstepper= false;
   
  private:
 
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void F01DetectorConstruction::ComputeCalorParameters()
{
     // Compute derived parameters of the calorimeter

     fZStartAbs = fZAbsorber-0.5*fAbsorberThickness;
     fZEndAbs   = fZAbsorber+0.5*fAbsorberThickness;

     G4cout << "-- Calorimeter Absorber z-coords:  Start= "
            << fZStartAbs << " End= " << fZEndAbs << G4endl;
}

#endif
