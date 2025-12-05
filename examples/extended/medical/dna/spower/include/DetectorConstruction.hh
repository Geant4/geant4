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
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"

class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction();
    virtual G4VPhysicalVolume* Construct();

    void SetRadius(G4double);
    void SetMaterial(G4String);
    void SetTrackingCut(G4double);

    inline G4double GetAbsorRadius() const { return fAbsorRadius; }

    inline G4Material* GetAbsorMaterial() const { return fAbsorMaterial; }

    inline G4double GetAbsorMass() const { return fLogicalAbsor->GetMass(); }

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
