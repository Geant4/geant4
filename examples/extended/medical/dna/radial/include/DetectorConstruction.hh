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
    void SetWorldOffsetLength(const G4double&);
    
    void SetCylinderLength(const G4double&);
    void SetCylinderThickness(const G4double&);
    void SetCylinderMinRadius(const G4double&);
    void SetCylinderMaxRadius(const G4double&);

    G4double GetWorldOffsetLength() const {return fWorldOffsetLength;}

    G4int GetCylinderNumber() const {return fCylinderNumber;}
    G4double GetCylinderThickness() const {return fCylinderThickness;}
    G4double GetCylinderLength() const {return fCylinderLength;}
    
  
  private:

    G4double fWorldRadius = 0;
    G4double fWorldOffsetLength = 0;

    G4double fCylinderLength = 0;
    G4double fCylinderThickness = 0;
    G4double fCylinderNumber = 0;
    G4double fCylinderMinRadius = 0;

    void DefineMaterials();

    G4LogicalVolume* fLogicWorld = nullptr;
    G4PVPlacement* fPhysiWorld = nullptr;

    G4Material* fWaterMaterial = nullptr;
    G4Material* fVacuumMaterial = nullptr;

    DetectorMessenger* fDetectorMessenger = nullptr;
};
#endif
