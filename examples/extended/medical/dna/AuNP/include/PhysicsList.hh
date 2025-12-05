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
/// \file PhysicsList.hh
/// \brief Definition of the PhysicsList class

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4ProcessManager.hh"
#include "G4VModularPhysicsList.hh"
#include <memory>

class DetectorConstruction;
class PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PhysicsList final : public G4VModularPhysicsList {
public:
  PhysicsList();
  ~PhysicsList() override = default;
  void SetPhysics4NP(const G4String &name);
  static void ConstructBosons();
  static void ConstructLeptons();
  static void ConstructBarions();
  static void ConstructGeneral() {};
  void ConstructEM();
  void ConstructParticle() override;
  void ConstructProcess() override;
  void SetCuts() override;
private:
  std::unique_ptr<G4VPhysicsConstructor> fEmDNAChemistryList;
  const DetectorConstruction *fpDetector = nullptr;
  std::unique_ptr<PhysicsListMessenger> fPhysMessenger;
  G4double fcutForGamma = 0.1 * CLHEP::nanometer;
  G4double fcutForElectron = 0.1 * CLHEP::nanometer;
  G4double fcutForPositron = 0.1 * CLHEP::nanometer;
  G4String fphysname;
};
#endif