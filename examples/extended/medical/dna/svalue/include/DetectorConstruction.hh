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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/include/DetectorConstruction.hh
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

  void SetTrackingCut(G4double);

  void SetCytoThickness(G4double);
  void SetCytoMaterial(G4String);
  
  void SetNuclRadius(G4double);
  void SetNuclMaterial(G4String);

  void SetWorldMaterial(G4String);

  virtual G4VPhysicalVolume* Construct();

  inline G4double GetCytoThickness() const
  {
    return fCytoThickness;
  }

  inline G4Material* GetCytoMaterial() const
  {
    return fCytoMaterial;
  }

  inline G4double GetCytoMass() const
  {
    return fLogicalCyto->GetMass();
  }

  inline G4double GetNuclRadius() const
  {
    return fNuclRadius;
  }

  inline G4Material* GetNuclMaterial() const
  {
    return fNuclMaterial;
  }

  inline G4double GetNuclMass() const
  {
    return fLogicalNucl->GetMass();
  }

  const G4LogicalVolume* GetNuclLogicalVolume() const
  {
    return fLogicalNucl;
  }

  const G4LogicalVolume* GetCytoLogicalVolume() const
  {
    return fLogicalCyto;
  }

  void PrintParameters() const;

private:

  void DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();

  G4double fTrackingCut;

  G4double fNuclRadius;
  G4Material* fNuclMaterial;
  
  G4double fCytoThickness;
  G4Material* fCytoMaterial;
  
  G4double fWorldRadius;
  G4Material* fWorldMaterial;

  G4VPhysicalVolume* fNucl;
  G4LogicalVolume* fLogicalNucl;

  G4VPhysicalVolume* fCyto;
  G4LogicalVolume* fLogicalCyto;
  
  G4VPhysicalVolume* fWorld;
  G4LogicalVolume* fLogicalWorld;

  DetectorMessenger* fDetectorMessenger;
};

#endif

