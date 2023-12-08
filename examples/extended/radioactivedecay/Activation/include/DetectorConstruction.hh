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

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction() override;

  public:

    G4VPhysicalVolume* Construct() override;

    G4Material* 
    MaterialWithSingleIsotope(G4String, G4String, G4double, G4int, G4int);

    void SetAbsorThickness(G4double);
    void SetAbsorSizeYZ   (G4double);
    void SetAbsorMaterial (G4String);

  public:  

   G4double GetAbsorThickness()    {return fAbsorThickness;};
   G4double GetAbsorSizeYZ()       {return fAbsorSizeYZ;};
   G4Material* GetAbsorMaterial()  {return fAbsorMaterial;};

   G4double GetWorldSizeX()   {return fWorldSizeX;};
   G4double GetWorldSizeYZ()  {return fWorldSizeYZ;};

   void PrintParameters();

  private:

   G4double           fAbsorThickness = 0.;
   G4double           fAbsorSizeYZ = 0.;
   G4Material*        fAbsorMaterial = nullptr;
   G4LogicalVolume*   fLAbsor = nullptr;

   G4double           fWorldSizeX = 0.;
   G4double           fWorldSizeYZ = 0.;
   G4Material*        fWorldMaterial = nullptr;
   G4VPhysicalVolume* fWorldVolume = nullptr;                        

   DetectorMessenger* fDetectorMessenger = nullptr;

  private:

   void               DefineMaterials();
   G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
