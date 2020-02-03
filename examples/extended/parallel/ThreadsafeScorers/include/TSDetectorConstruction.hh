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
/// \file parallel/ThreadsafeScorers/include/TSDetectorConstruction.hh
/// \brief Definition of the TSDetectorConstruction class
//
//
//
//
/// Construction of a target material (default = boron) surrounded by a
///     casing material (default = water) and a vacuum world (default =
///     target and casing fill world). The target + casing is brick
///     geometry with fTargetSections defining the number of divisions
///     in each dimension. The end sections in each dimension
///     is set to the casing. So a fTargetSections = G4ThreeVector(3, 3, 3)
///     would be one section of boron and 8 sections of water.
/// The idea behind this geometry is just to create a simple geometry that
///     scatters and produces a lot neutrons with a minimal number of sections
///     (i.e. coarse meshing) such that the contention in operating on
///     the atomic hits maps is higher and round-off errors in the
///     thread-local hits maps are detectable (printed out in TSRunAction)
///     from the sheer number of floating point sum operations.
/// Two scorers are implemented: EnergyDeposit and Number of steps
///     The energy deposit is to (possibly) show the round-off error seen
///     with thread-local hits maps. The # of steps scorer is to verify
///     the thread-safe and thread-local hits maps provide the same results.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef tsdetectorconstruction_hh
#define tsdetectorconstruction_hh 1


#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"

#include <map>
#include <set>

class G4Box;
class G4Tubs;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class TSDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    typedef std::map<G4String, G4Material*> MaterialCollection_t;
    typedef std::set<G4LogicalVolume*> ScoringVolumes_t;

public:
    TSDetectorConstruction();
    virtual ~TSDetectorConstruction();

    static TSDetectorConstruction* Instance();

public:
    G4VPhysicalVolume* Construct();
    inline const G4ThreeVector& GetWorldDimensions() const { return fWorldDim; }
    inline const ScoringVolumes_t& GetScoringVolumes() const
    { return fScoringVolumes; }
    inline const G4String& GetMFDName() const { return fMfdName; }
    inline G4int GetTotalTargets() const 
    { return fTargetSections.x() * fTargetSections.y() * fTargetSections.z(); }

protected:
    virtual MaterialCollection_t ConstructMaterials();
    virtual G4VPhysicalVolume* ConstructWorld(const MaterialCollection_t&);
    virtual void ConstructSDandField();

private:
    static TSDetectorConstruction* fgInstance;
    G4VPhysicalVolume* fWorldPhys;
    ScoringVolumes_t fScoringVolumes;
    G4String fWorldMaterialName;
    G4String fTargetMaterialName;
    G4String fCasingMaterialName;
    G4ThreeVector fWorldDim;
    G4ThreeVector fTargetDim;
    G4ThreeVector fTargetSections;
    G4String fMfdName;
};

#endif
