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

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include <memory>
class G4VPhysicalVolume;
class G4Material;
class G4MultiFunctionalDetector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction final : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction() override = default;
  void SetNReplicaR(G4int);
  void SetNReplicaAzm(G4int);
  void SetAbsRadius(G4double);
  void SetAbsMaterial(const G4String&);
  void SetNPRadius(G4double);
  void SetNPMaterial(const G4String&);
  void SetTrackingCut(G4double);
  G4VPhysicalVolume *Construct() override;
  G4double GetNReplicaR() const { return fNreplicaR; }
  G4double GetNReplicaAzm() const { return fNreplicaAzm; }
  G4double GetAbsRadius() const { return fAbsRadius; }
  G4double GetNPRadius() const { return fNPRadius; }
  G4Material *GetNPMaterial() const { return fNPMaterial; }
  G4double GetTrackingCut() const { return fTrackingCut; }
  G4double GetNPMass() const { return fLogicalNP->GetMass(); }
  G4MultiFunctionalDetector *GetMFDetector() const { return fMFD; }
  void PrintParameters() const;
  G4Region *GetTargetRegion() const { return fRegion;}
private:
  static void DefineMaterials();
  G4VPhysicalVolume *ConstructVolumes();
  G4int fNreplicaR = 1000;
  G4int fNreplicaAzm = 360;
  G4double fTrackingCut = 10.0 * CLHEP::eV; // default tracking cut
  G4double fNPRadius = 50 * CLHEP::nm;
  G4double fAbsRadius = 100000 * CLHEP::nm + fNPRadius;
  G4Material *fNPMaterial = nullptr;
  G4Material *fAbsMaterial = nullptr;
  G4VPhysicalVolume *pWorld = nullptr;
  G4VPhysicalVolume *fNP = nullptr;
  G4VPhysicalVolume *fAbs = nullptr;
  G4LogicalVolume *fLogicalNP = nullptr;
  G4LogicalVolume *fLogicalAbs = nullptr;
  std::unique_ptr<DetectorMessenger> fDetectorMessenger;
  G4MultiFunctionalDetector *fMFD = nullptr;
  G4Region *fRegion = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
