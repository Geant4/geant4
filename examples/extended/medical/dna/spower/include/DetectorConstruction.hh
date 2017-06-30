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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"

class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction: public G4VUserDetectorConstruction
{

public:

  DetectorConstruction();
  ~DetectorConstruction();

  void SetRadius(G4double);
  
  void SetMaterial(G4String);
  
  void SetTrackingCut(G4double);
  
  virtual G4VPhysicalVolume* Construct();

  inline G4double GetAbsorRadius() const
  {
    return fAbsorRadius;
  }
  
  inline G4Material* GetAbsorMaterial() const
  {
    return fAbsorMaterial;
  }
  
  inline G4double GetAbsorMass() const
  {
    return fLogicalAbsor->GetMass();
  }
  
  void PrintParameters() const;

private:

  void DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();

  G4double fAbsorRadius;
  G4Material* fAbsorMaterial;
  G4VPhysicalVolume* fAbsor;
  G4LogicalVolume* fLogicalAbsor;
  DetectorMessenger* fDetectorMessenger;
};

#endif

