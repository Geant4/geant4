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
/// \file electromagnetic/TestEm5/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4Box;
class G4VPhysicalVolume;
class G4Material;
class G4MaterialCutsCouple;
class G4UniformMagField;
class DetectorMessenger;
class G4GlobalMagFieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  explicit DetectorConstruction();
          ~DetectorConstruction() override;

  void SetAbsorberMaterial (const G4String&);
  void SetAbsorberThickness(G4double);
  void SetAbsorberSizeYZ   (G4double);

  void SetAbsorberXpos(G4double);

  void SetWorldMaterial(const G4String&);
  void SetWorldSizeX   (G4double);
  void SetWorldSizeYZ  (G4double);

  void SetMagField(G4double);

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  void PrintGeomParameters();

  const G4Material* GetAbsorberMaterial() const {return fAbsorberMaterial;};
  G4double GetAbsorberThickness() const         {return fAbsorberThickness;};
  G4double GetAbsorberSizeYZ() const            {return fAbsorberSizeYZ;};

  G4double GetAbsorberXpos() const              {return fXposAbs;};
  G4double GetxstartAbs() const                 {return fXstartAbs;};
  G4double GetxendAbs() const                   {return fXendAbs;};

  const G4Material* GetWorldMaterial() const    {return fWorldMaterial;};
  G4double GetWorldSizeX() const                {return fWorldSizeX;};

  const G4VPhysicalVolume* GetAbsorber() const  {return fPhysiAbsorber;};

private:

  void DefineMaterials();
  void ComputeGeomParameters();
  void ChangeGeometry();

  G4Material*        fAbsorberMaterial = nullptr;
  G4double           fAbsorberThickness = 0.;
  G4double           fAbsorberSizeYZ = 0.;

  G4double           fXposAbs = 0.;
  G4double           fXstartAbs = 0., fXendAbs = 0.;

  G4Material*        fWorldMaterial = nullptr;
  G4double           fWorldSizeX = 0.;
  G4double           fWorldSizeYZ = 0.;

  G4Box*             fSolidWorld = nullptr;
  G4LogicalVolume*   fLogicWorld = nullptr;
  G4VPhysicalVolume* fPhysiWorld = nullptr;

  G4Box*             fSolidAbsorber = nullptr;
  G4LogicalVolume*   fLogicAbsorber = nullptr;
  G4VPhysicalVolume* fPhysiAbsorber = nullptr;
     
  DetectorMessenger* fDetectorMessenger = nullptr;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

