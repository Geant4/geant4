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
/// \file electromagnetic/TestEm2/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4Tubs;
class G4LogicalVolume;
class DetectorMessenger;
class G4GlobalMagFieldMessenger;

const G4int kMaxBin = 500;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
 ~DetectorConstruction() override;

public:

  void SetMaterial(const G4String&);
  void SetLBining (G4ThreeVector);
  void SetRBining (G4ThreeVector);

  G4VPhysicalVolume* Construct() override;

  void ConstructSDandField() override;

  const
  G4VPhysicalVolume* GetEcal() const    {return fPhysiEcal;};
  const G4Material* GetMaterial() const {return fMaterial;};

  // Subdivision of absorber
  G4int    GetnLtot() const      {return fNLtot;};
  G4int    GetnRtot() const      {return fNRtot;};
  G4double GetdLradl() const     {return fDLradl;};
  G4double GetdRradl() const     {return fDRradl;};
  G4double GetdLlength() const   {return fDLlength;};
  G4double GetdRlength() const   {return fDRlength;};     
  G4double GetfullLength() const {return fEcalLength;};
  G4double GetfullRadius() const {return fEcalRadius;};

private:

  void DefineMaterials();
  void UpdateParameters();

  G4int    fNLtot = 40,    fNRtot = 50;     // nb of bins: longitudinal and radial
  G4double fDLradl = 0.5,  fDRradl = 0.1;   // bin thickness (in radl unit)
  G4double fDLlength = 0., fDRlength = 0.;  // bin thickness (in length unit)

  G4Material* fMaterial = nullptr;       //pointer to the material
    
  G4double fEcalLength = 0.;             //full length of the Calorimeter
  G4double fEcalRadius = 0.;             //radius  of the Calorimeter

  G4Tubs*            fSolidEcal = nullptr;  //pointer to the solid calorimeter
  G4LogicalVolume*   fLogicEcal = nullptr;  //pointer to the logical calorimeter
  G4VPhysicalVolume* fPhysiEcal = nullptr;  //pointer to the physical calorimeter

  DetectorMessenger* fDetectorMessenger = nullptr;  //pointer to the Messenger

  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger = nullptr;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

