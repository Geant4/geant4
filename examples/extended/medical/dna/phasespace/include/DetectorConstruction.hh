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

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"

class DetectorMessenger;
class G4LogicalVolume;
class G4PVPlacement;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;
    
    G4VPhysicalVolume* Construct() override;
    
    void SetMaterial(const G4String&);
    void SetWorldSize(G4double); 
    void SetScorerRadius(G4double); 

  public:
    
    G4Material* MaterialWithDensity(G4String, G4double); 

    G4double GetWorldSize() {return fWorldSize;};
    G4double GetScorerRadius() {return fScorerRadius;};
    const G4LogicalVolume* GetLogicalWorld() const {return fpLogicWorld;};
    const G4LogicalVolume* GetLogicalScorer() const {return fpLogicScorer;}; 
  
  private:
   
    G4double fWorldSize = 0.;
    G4double fScorerRadius = 0.;
    
    void DefineMaterials();

    DetectorMessenger* fpDetectorMessenger = nullptr;
    G4Material* fpWaterMaterial = nullptr;
    G4LogicalVolume* fpLogicWorld = nullptr;
    G4PVPlacement* fpPhysiWorld = nullptr;

    G4LogicalVolume* fpLogicScorer = nullptr;
    G4PVPlacement* fpPhysiScorer = nullptr;
};

#endif
