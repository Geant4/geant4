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
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <CLHEP/Units/SystemOfUnits.h>
#include "G4VUserDetectorConstruction.hh"

#include <memory>

class G4Material;
class DetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction()override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;
    
    void SetTrackingCut(const G4double&);
    void SetMaxStepSize(const G4double&);
    void PrintParameters() const;
    void SetTrackerSDRadius(const G4double&);

  private:
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    std::unique_ptr<DetectorMessenger> fDetectorMessenger;
    G4Material* fpWaterMaterial = nullptr;
    G4double fpTrackingCut = 11.* CLHEP::eV;
    G4double fpMaxStepSize = DBL_MAX;
    G4double fRadius = 0.;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
