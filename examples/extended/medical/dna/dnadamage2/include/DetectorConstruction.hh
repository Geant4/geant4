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
// DNADAMAGE2 example is derived from the chem6 example
// chem6 example authors: W. G. Shin and S. Incerti (CENBG, France)
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DNADAMAGE2_DetectorConstruction_h
#define DNADAMAGE2_DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4UImessenger.hh"
#include "G4Orb.hh"
#include "G4MoleculeGun.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "PhysGeoImport.hh"
#include "StackingAction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction,
                             public G4UImessenger
{
public:
  DetectorConstruction();
  ~DetectorConstruction() override;
  void SetNewValue(G4UIcommand*,G4String) override;
  void SetSize(G4double);

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;
  void ReadOffsetFile(G4String);
  void AddDNAInformation(G4int, G4ThreeVector);
  void SetStacking(StackingAction* stack) {fpStacking = stack;}
  StackingAction* GetStacking() {return fpStacking;}

  std::vector<G4String> GetDNANames() {return fDNANames;}
  std::vector<G4ThreeVector> GetDNAPositions() {return fDNAPositions;}
  std::vector<std::vector<G4int>> GetDNADetails() {return fDNADetails;}

private:
  G4UIdirectory* fDetDir = nullptr;
  G4UIcmdWithAString* fpOffSetFileUI = nullptr;
  G4UIcmdWithAString* fpPlasmidFile  = nullptr;
  G4UIcmdWithAnInteger* fpPlasmidNbUI = nullptr;
  G4UIcmdWithADoubleAndUnit* fSizeCmd = nullptr;
  G4UIcmdWithABool* fpUseDNA = nullptr;

  G4Orb* fPlasmidEnvelope = nullptr;
  G4double fWorldSize = 1 * um;

  G4int fNbOfPlasmids = 0;
  G4String fPlasmidFile = "VoxelStraight.fab2g4dna";
  G4bool fUseDNAVolumes = false;
  std::vector<G4ThreeVector> fVOffset;

  std::vector<G4String> fDNANames;
  std::vector<G4ThreeVector> fDNAPositions;
  std::vector<std::vector<G4int>> fDNADetails;

  std::vector<G4String> fSampleDNANames;
  std::vector<G4ThreeVector> fSampleDNAPositions;
  std::vector<std::vector<G4int>> fSampleDNADetails;

  StackingAction* fpStacking = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
