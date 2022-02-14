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
//

#pragma once
#include <map>
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include <memory>
#include <vector>
#include <unordered_map>
#include "DNAVolumeType.hh"
#include "G4MoleculeGun.hh"

class G4MoleculeGun;
class G4VPhysicalVolume;
class G4VSolid;
class G4Material;
class G4LogicalVolume;
class G4Orb;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAParser
{
public:
    using GeoData = std::unordered_map<const G4VPhysicalVolume*, 
                                       DNAVolumeType>;
    DNAParser();
    ~DNAParser();
    G4LogicalVolume* CreateLogicVolume();
    G4double GetVoxelSize();      
    
    std::unique_ptr<G4MoleculeGun> ReleaseMoleculeGun()
    {
        return std::move(fpGun);
    }
    
    GeoData ReleaseGeoData()
    {
        return std::move(fGeometryMap);
    }
    
    struct Molecule;
    
    void ParseFile(const std::string&);
private:
    G4double fSize;
    std::string fGeoName;
    std::map<G4String, G4double> fRadiusMap;
    std::map<G4String, G4double> fWaterRadiusMap;
    std::vector<Molecule> fMolecules;
    G4Material* fpWater;
    std::unique_ptr<G4MoleculeGun> fpGun;
    std::map<std::string, DNAVolumeType> fEnumMap;
    GeoData fGeometryMap;
    void EnumParser();
    G4VSolid* CreateCutSolid(G4Orb*,
                             Molecule&,
                             std::vector<Molecule>&,
                             G4bool);
};
