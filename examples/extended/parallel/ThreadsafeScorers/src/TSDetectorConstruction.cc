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
/// \file parallel/ThreadsafeScorers/src/TSDetectorConstruction.cc
/// \brief Implementation of the TSDetectorConstruction class
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




#include "TSDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSNofStep.hh"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSDetectorConstruction* TSDetectorConstruction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSDetectorConstruction* TSDetectorConstruction::Instance()
{
  return fgInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSDetectorConstruction::TSDetectorConstruction()
: fWorldPhys(0),
  fWorldMaterialName("G4_Galactic"),
  fTargetMaterialName("G4_B"),
  fCasingMaterialName("G4_WATER"),
  fWorldDim(G4ThreeVector(0.5*m, 0.5*m, 0.5*m)),
  fTargetDim(G4ThreeVector(0.5*m, 0.5*m, 0.5*m)),
  fTargetSections(G4ThreeVector(5, 5, 5)),
  fMfdName("Target_MFD")
{
  fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSDetectorConstruction::~TSDetectorConstruction()
{
  fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* TSDetectorConstruction::Construct()
{
    return ConstructWorld(ConstructMaterials());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSDetectorConstruction::MaterialCollection_t
TSDetectorConstruction::ConstructMaterials()
{
    MaterialCollection_t materials;
    G4NistManager* nist = G4NistManager::Instance();

    materials["World"] = nist->FindOrBuildMaterial(fWorldMaterialName);
    materials["Target"] = nist->FindOrBuildMaterial(fTargetMaterialName);
    materials["Casing"] = nist->FindOrBuildMaterial(fCasingMaterialName);

    return materials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume*
TSDetectorConstruction::ConstructWorld(const MaterialCollection_t& materials)
{
    G4UserLimits* steplimit = new G4UserLimits(0.1*(fTargetDim.z()
                                                    /fTargetSections.z()));
    G4bool check_overlap = false;

    G4Box* world_solid = new G4Box("World",
                                   0.5*fWorldDim.x(),
                                   0.5*fWorldDim.y(),
                                   0.5*fWorldDim.z());
    G4LogicalVolume* world_log = new G4LogicalVolume(world_solid,
                                                     materials.find("World")
                                                     ->second,
                                                     "World");
    fWorldPhys = new G4PVPlacement(0,
                                   G4ThreeVector(0.),
                                   "World",
                                   world_log,
                                   0, false, 0, check_overlap);

    G4int nz = fTargetSections.z();
    G4int ny = fTargetSections.y();
    G4int nx = fTargetSections.x();

    // spacing between sections
    G4double sx = fTargetDim.x()/fTargetSections.x();
    G4double sy = fTargetDim.y()/fTargetSections.y();
    G4double sz = fTargetDim.z()/fTargetSections.z();

    //G4cout << "World has dimensions : "
    //<< G4BestUnit(fWorldDim, "Length") << G4endl;

    //------------------------------------------------------------------------//
    // Set Visual Attributes
    //------------------------------------------------------------------------//
    G4VisAttributes* red   = new G4VisAttributes(G4Color(1., 0., 0., 1.0));
    G4VisAttributes* green = new G4VisAttributes(G4Color(0., 1., 0., 0.25));
    G4VisAttributes* blue  = new G4VisAttributes(G4Color(0., 0., 1., 0.1));
    G4VisAttributes* white = new G4VisAttributes(G4Color(1., 1., 1., 1.));

    white->SetVisibility(true);
    red->SetVisibility(true);
    green->SetVisibility(true);
    blue->SetVisibility(true);

    white->SetForceWireframe(true);
    red->SetForceSolid(true);
    green->SetForceSolid(true);
    blue->SetForceSolid(true);

    world_log->SetVisAttributes(white);

    for(G4int k = 0; k < nz; ++k)
        for(G4int j = 0; j < ny; ++j)
            for(G4int i = 0; i < nx; ++i)
            {
                // displacement of section
                G4double dx
                = 0.5*sx + static_cast<G4double>(i)*sx - 0.5*fWorldDim.x();
                G4double dy
                = 0.5*sy + static_cast<G4double>(j)*sy - 0.5*fWorldDim.y();
                G4double dz
                = 0.5*sz + static_cast<G4double>(k)*sz - 0.5*fWorldDim.z();
                G4ThreeVector td = G4ThreeVector(dx, dy, -dz);
                // make unique name
                std::stringstream ss_name;
                ss_name << "Target_" << i << "_" << j << "_" << k;

                G4Box* target_solid = new G4Box(ss_name.str(),
                                                0.5*sx,
                                                0.5*sy,
                                                0.5*sz);


                G4Material* target_material = 0;
                G4bool is_casing = true;

                if(j == 0 || j+1 == ny || i == 0 || i+1 == nx ||
                   (nz > 1 && (k == 0 || k+1 == nz)))
                    target_material = materials.find("Casing")->second;
                else {
                    target_material = materials.find("Target")->second;
                    is_casing = false;
                }

                G4LogicalVolume* target_log =
                new G4LogicalVolume(target_solid,
                                    target_material,
                                    ss_name.str());

                target_log->SetUserLimits(steplimit);

                new G4PVPlacement(0,
                                  td,
                                  ss_name.str(),
                                  target_log,
                                  fWorldPhys,
                                  true,
                                  k*nx*ny + j*nx + i,
                                  check_overlap);

                fScoringVolumes.insert(target_log);

                if(is_casing)
                    target_log->SetVisAttributes(blue);
                else
                {
                    // making a checkerboard for kicks...
                    G4bool even_z = (k%2 == 0) ? true : false;
                    G4bool even_y = (j%2 == 0) ? true : false;
                    G4bool even_x = (i%2 == 0) ? true : false;

                    G4VisAttributes* theColor = nullptr;

                    if((even_z))
                    {
                      if((even_y && even_x) || (!even_y && !even_x))
                        theColor = red;
                      else
                        theColor = green;
                    } else // ! even_z
                    {
                      if((!even_y && even_x) || (even_y && !even_x))
                        theColor = red;
                      else
                        theColor = green;
                    }

                    target_log->SetVisAttributes(theColor);
                }

            }



    return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TSDetectorConstruction::ConstructSDandField()
{
    //------------------------------------------------//
    //            MultiFunctionalDetector             //
    //------------------------------------------------//
    // Define MultiFunctionalDetector with name.
    G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector(fMfdName);
    G4SDManager::GetSDMpointer()->AddNewDetector(MFDet);
    G4VPrimitiveScorer* edep = new G4PSEnergyDeposit("EnergyDeposit");
    MFDet->RegisterPrimitive(edep);
    G4VPrimitiveScorer* nstep = new G4PSNofStep("NumberOfSteps");
    MFDet->RegisterPrimitive(nstep);

    // add scoring volumes
    for(auto ite : fScoringVolumes)
    {
        SetSensitiveDetector(ite, MFDet);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


