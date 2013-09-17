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
// and the DNA geometry given in the Geom_DNA example 
// shall cite the following Geant4-DNA collaboration publications:
// [1] NIM B 298 (2013) 47-54
// [2] Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4Orb.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4UnionSolid.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();

  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();
                         
private:
   
  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;
  G4Box*             solidWorld;

  G4VPhysicalVolume* physiTin;
  G4LogicalVolume*   logicTin;
  G4Box*             solidTin;

  G4VPhysicalVolume* physiNucleus;
  G4LogicalVolume*   logicNucleus;
  G4Ellipsoid*       solidNucleus;

  G4Tubs*            solidBoxros;
  G4VPhysicalVolume* physiBox[48];
  G4VPhysicalVolume* physiBoxros[48];
  G4LogicalVolume*   logicBoxros;

  G4RotationMatrix*  rotTrap;
  G4ThreeVector      posTrap;

  G4VPhysicalVolume* physiEnv;
  G4LogicalVolume*   logicEnv;

  G4LogicalVolume*   logicHistone;
  G4Tubs*            solidHistone;

  G4VPhysicalVolume* physiBp1[200];
  G4LogicalVolume*   logicBp1;
  G4Orb*             solidBp1;

  G4VPhysicalVolume* physiBp2[200];
  G4LogicalVolume*   logicBp2;
  G4Orb*             solidBp2;

  G4Orb*             solidSugar1;
  G4Orb*             solidSugar2;
  G4Orb*             solidSugar3;
  G4Orb*             solidSugar4;

  G4LogicalVolume*   logicSphere3;
  G4Orb*             solidSphere3;

  G4LogicalVolume*   logicSphere4;
  G4Orb*             solidSphere4;

  G4UnionSolid* uniDNA;
  G4UnionSolid* uniDNA2;


  G4Material*        waterMaterial;


  void DefineMaterials();
  G4VPhysicalVolume* ConstructDetector();     
  void LoadChromosome(const char* filename, G4VPhysicalVolume* physiBox);
};
#endif
