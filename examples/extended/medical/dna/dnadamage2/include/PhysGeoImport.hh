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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Update: H. Tran (IRSN, France) :20/12/2018
//         J. Naoki D. Kondo (UCSF, US): 10/10/2021
//
/// \file PhysGeoImport.hh
/// \brief Definition of the plasmid load methods for the geometry

#ifndef DNADAMAGE2_GeoImport_h
#define DNADAMAGE2_GeoImport_h 1

#include <map>
#include <fstream>
#include <algorithm>

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Orb.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include <memory>
#include "G4H2O.hh"
#include "G4Electron_aq.hh"

class G4VPhysicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct Molecule
{
    Molecule(G4String name, G4int copyNumber, G4ThreeVector position, 
        G4double radius, G4double waterRadius, std::string material, G4int strand)
    {
        fName = name;
        fMaterial = material;
        fCopyNumber = copyNumber;
        fPosition = position;
        fRadius = radius;
        fRadiusWater = waterRadius;
        fStrand = strand;
    }

    G4String fName     = "none";
    G4String fMaterial = "none";

    G4int fCopyNumber = -1;
    G4int fStrand     = -1;

    G4ThreeVector fPosition = G4ThreeVector();

    G4double fRadius = 1;
    G4double fRadiusWater = 1;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysGeoImport
{
public:
    PhysGeoImport();
    ~PhysGeoImport() = default;

    G4LogicalVolume* CreateLogicVolumeXYZ(G4String fileName);

    std::vector<G4String> GetMoleculesNames() {return fSampleDNANames;}
    std::vector<G4ThreeVector> GetMoleculesPositions() {return fSampleDNAPositions;}
    std::vector<std::vector<G4int>> GetMoleculesDetails() {return fSampleDNADetails;}

private:
    std::string fGeoName = "VoxelStraight";

    std::map<G4String, G4double> fRadiusMap;
    std::map<G4String, G4double> fWaterRadiusMap;

    std::vector<Molecule> fMolecules;

    // Materials
    G4Material* fpWater = nullptr;
    G4Material* fEnvelopeWater = nullptr;

    void ReadFile(G4String fileName);

    G4double fOffsetX = 0;
    G4double fOffsetY = 0;
    G4double fOffsetZ = 0;
    G4double fXMin =  1000;
    G4double fXMax = -1000;
    G4double fYMin =  1000;
    G4double fYMax = -1000;
    G4double fZMin =  1000;
    G4double fZMax = -1000;
    std::ofstream fOutDNA;

    std::vector<G4ThreeVector> fVertexes;

    std::vector<G4String> fSampleDNANames;
    std::vector<G4ThreeVector> fSampleDNAPositions;
    std::vector<std::vector<G4int>> fSampleDNADetails;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
