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
// Med. Phys. 51 (2024) 5873-5889
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4VUserDetectorConstruction.hh"

class DetectorMessenger;
class PhysicsList;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction(PhysicsList*);
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;

    void SetMaterial(const G4String&);
    void SetWorldRadius(const G4double&);
    void SetWorldLength(const G4double&);
    void SetThicknessCylinders(const G4double&);
    void SetMinRadiusCylinders(const G4double&);
    void SetMaxRadiusCylinders(const G4double&);

    G4int GetNumberCylinders() const {return fNumberCylinders;}
    G4double GetThicknessCylinders() const {return fThicknessCylinders;}
    G4double GetWorldLength() const {return fWorldLength;}

  private:

    G4double fWorldRadius = 0;
    G4double fWorldLength = 0;

    G4double fThicknessCylinders = 0;
    G4double fNumberCylinders = 0;
    G4double fMinRadiusCylinders = 0;

    void DefineMaterials();

    G4LogicalVolume* fLogicWorld = nullptr;
    G4PVPlacement* fPhysiWorld = nullptr;

    G4Material* fWaterMaterial = nullptr;
    G4Material* fVacuumMaterial = nullptr;

    DetectorMessenger* fDetectorMessenger = nullptr;
};
#endif
