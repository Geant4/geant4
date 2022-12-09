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
// -------------------------------------------------------------
//  =============== Begin Documentation Comments ===============
//!
//! \file       FFDetectorConstruction.cc
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Implementation of the FFDetectorConstruction class
//!
//! \details    The model simulated is based off a subcritical assembly design
//!             with 20% enriched meat
//!
//  ================ End Documentation Comments ================
//
//  Modified:
//
//  23-06-14                                              BWendt
//  Fixed issue with the automatic placement of the meat not working. Solution
//      was to use the correct units "inch" in the y-direction as well.
//  Coincidentally eliminated the need to using the 'std::abs()" from the
//      "cmath" library.
//  Implemented method "PlaceFuelPlates"
//
// -------------------------------------------------------------

#include "globals.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

#include "FFDetectorConstruction.hh"


static const G4double inch = 2.54 * cm;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
FFDetectorConstruction::
FFDetectorConstruction()
:   G4VUserDetectorConstruction()
{
    DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* FFDetectorConstruction::
Construct()
{
    G4ThreeVector position;
#ifdef NDEBUG
    G4bool const overlapChecking = false;
#else
    G4bool const overlapChecking = true;
#endif // NDEBUG
    
    //
    // Create the world
    //
    const G4double worldSize = 40.0 * inch;
    G4Box* const solidWorld = new G4Box("World",        // the name
                                        worldSize,      // x size
                                        worldSize,      // y size
                                        worldSize);     // z size
    G4LogicalVolume* const logicalWorld
        = new G4LogicalVolume(solidWorld,               // the solid volume
                              fAir,                     // the material
                              solidWorld->GetName());   // the name
    // Center at the origin
    position.set(0.0, 0.0, 0.0); 
    G4VPhysicalVolume* const physicalWorld
        = new G4PVPlacement(NULL,                       // no rotation
                            position,                   // must be at origin
                            logicalWorld,               // the logical volume
                            logicalWorld->GetName(),    // the name
                            NULL,                       // no mother volume
                            false,                      // no boolean ops
                            0,                          // copy number
                            overlapChecking);           // check for overlaps
    
    //
    // Create the graphite pile that the subcritical assembly rests on.
    //
    const G4double floorH = 30.0 * inch;
    const G4ThreeVector floorPosition(0.0, 0.0, 0.0);
    G4Box* const solidFloor = new G4Box("Floor",        // the name
                                        worldSize,      // x size
                                        worldSize,      // y size
                                        floorH * 0.5);  // z size
    G4LogicalVolume* const logicalFloor
        = new G4LogicalVolume(solidFloor,               // the solid volume
                              fGraphite,                // the material
                              solidFloor->GetName());   // the name
    // Shift down so the top is at the origin
    position.set(0.0, 0.0, -floorH * 0.5); 
    new G4PVPlacement(NULL,                             // no rotation
                      position,                         // position
                      logicalFloor,                     // the logical volume
                      logicalFloor->GetName(),          // the name
                      logicalWorld,                     // the mother volume
                      false,                            // no boolean ops
                      0,                                // copy number
                      overlapChecking);                 // check for overlaps
    
    //
    // Create the tank
    //
    const G4double tankWallThickness = 0.25 * inch;
    const G4double tankOR = 18.0 * inch;
    const G4double tankH = 39.0 * inch;
    G4Tubs* const solidTank
        = new G4Tubs("Tank_Wall",                       // the name
                     0.0,                               // inner radius
                     tankOR,                            // outer radius
                     tankH * 0.5,                       // half height
                     0.0 * deg,                         // start angle
                     360.0 * deg);                      // end angle
    G4LogicalVolume* const logicalTank
        = new G4LogicalVolume(solidTank,                // the solid volume
                              fAluminum,                // the material
                              solidTank->GetName());    // the name
    // Shift up so the base is at the origin
    position.set(0.0, 0.0, tankH * 0.5);
    new G4PVPlacement(NULL,                             // no rotation
                      position,                         // shift up
                      logicalTank,                      // the logical volume
                      logicalTank->GetName(),           // the name
                      logicalWorld,                     // the mother volume
                      false,                            // no boolean ops
                      0,                                // copy number
                      overlapChecking);                 // check for overlaps
    // Top 3 inches are air
    const G4double tankAirH = 3.0 * inch;
    G4Tubs* const solidTankAir
        = new G4Tubs("Tank_Air",                        // the name
                     0.0,                               // inner radius
                     tankOR - tankWallThickness,        // outer radius
                     tankAirH * 0.5,                    // half height
                     0.0 * deg,                         // start angle
                     360.0 * deg);                      // end angle
    G4LogicalVolume* const logicalTankAir
        = new G4LogicalVolume(solidTankAir,             // the solid volume
                              fAir,                     // the material
                              solidTankAir->GetName()); // the name
    // Shift up so that the top of the air is the same as the top of the tank
    position.set(0.0, 0.0, (tankH - tankAirH) * 0.5);
    new G4PVPlacement(NULL,                             // no rotation
                      position,                         // shift ip
                      logicalTankAir,                   // the logical volume
                      logicalTankAir->GetName(),        // the name
                      logicalTank,                      // the mother volume
                      false,                            // no boolean ops
                      0,                                // copy number
                      overlapChecking);                 // check for overlaps
    // Fill remaining area with water
    const G4double tankH2OH = (tankH - (tankAirH + tankWallThickness));
    G4Tubs* const solidTankH2O
        = new G4Tubs("Tank_H2O",                        // the name
                     0.0,                               // inner radius
                     tankOR - tankWallThickness,        // outer radius
                     tankH2OH * 0.5,                    // half height
                     0.0 * deg,                         // start angle
                     360.0 * deg);                      // end angle
    G4LogicalVolume* const logicalTankH2O
        = new G4LogicalVolume(solidTankH2O,             // the solid volume
                              fAluminum,                // the material
                              solidTankH2O->GetName()); // the name
    // Shift up so that the top of the water is at the bottom of the air
    const G4double centerOfH2O = (tankH - tankH2OH) * 0.5 - tankAirH;
    position.set(0.0, 0.0, centerOfH2O);
    new G4PVPlacement(NULL,                             // no rotation
                      position,                         // shift to origin
                      logicalTankH2O,                   // the logical volume
                      logicalTankH2O->GetName(),        // the name
                      logicalTank,                      // the mother volume
                      false,                            // no boolean ops
                      0,                                // copy number
                      overlapChecking);                 // check for overlaps
    
    //
    // Fuel plates
    //
    const G4double plateX = 3.0 * inch;
    const G4double plateY = 0.08 * inch;
    const G4double plateZ = 26.0 * inch;
    const G4double meatX = 2.75 * inch;
    const G4double meatY = 0.04 * inch;
    const G4double meatZ = 24.0 * inch;
    const G4double xSpacing = 5.0 * inch;
    const G4double ySpacing = 0.3 * inch;
    const G4double plateRadius = 12.0 * inch;
    // Define the aluminim claddiing
    G4Box* const solidPlate
        = new G4Box("Plate_Cladding",                   // the name
                    plateX * 0.5,                       // x size
                    plateY * 0.5,                       // y size
                    plateZ * 0.5);                      // z size
    G4LogicalVolume* const logicalPlate
        = new G4LogicalVolume(solidPlate,               // the solid volume
                              fAluminum,                // the material
                              solidPlate->GetName());   // the name
    // Place the meat inside the cladding
    G4Box* const solidMeat
        = new G4Box("Plate_Meat",                       // the name
                    meatX * 0.5,                        // x size
                    meatY * 0.5,                        // y size
                    meatZ * 0.5);                       // z size
    G4LogicalVolume* const logicalMeat
        = new G4LogicalVolume(solidMeat,                // the solid volume
                              fUO2_20E,                 // the material
                              solidMeat->GetName());    // the name
    // The meat goes into the exact center of the plate
    position.set(0.0, 0.0, 0.0); 
    new G4PVPlacement(NULL,                             // no rotation
                      position,                         // position
                      logicalMeat,                      // the logical volume
                      logicalMeat->GetName(),           // the name
                      logicalPlate,                     // the mother volume
                      false,                            // no boolean ops
                      0,                                // copy number
                      overlapChecking);                 // check for overlaps
    // The plate will be centered in the z-direction within the water
    // Simulate a subcritical assembly loading within a radius of 12 inches
    bool placeMe;
    
    position.setZ(0.0);
    fCopyNumber = 0;
    for(double x = 0.0;
        x <= plateRadius;
        x += xSpacing)
    {
        // 5 rows of plates
        for(double y = 0.0;
            y <= plateRadius;
            y += ySpacing)
        {
            placeMe = false;
            
            // Fuel plate must be completely within the radius to be placed
            if(std::sqrt(x * x + y * y) < plateRadius)
            {
                // Leave a 1 inch radius opening in the middle for the neutron
                // source
                if(std::sqrt(x * x + y * y) > 1.0 * inch)
                {
                    placeMe = true;
                }
            }
            
            if(placeMe)
            {
                PlaceFuelPlate(x,
                               y,
                               logicalPlate,
                               logicalTankH2O);
                PlaceFuelPlate(x,
                               -y,
                               logicalPlate,
                               logicalTankH2O);
                if(x > 0.0)
                {
                    PlaceFuelPlate(-x,
                                   y,
                                   logicalPlate,
                                   logicalTankH2O);
                    PlaceFuelPlate(-x,
                                   -y,
                                   logicalPlate,
                                   logicalTankH2O);
                }
            }
        }
    }
    G4cout << fCopyNumber << " plates were added to the subcritical assembly"
           << G4endl;
    
    //
    // Neutron Source
    //
    // TODO create the AmBe material in DefineMaterials() and use it here
    //      For now steel is used, but the logical volume is used in the
    //      PrimaryGeneratorAction to know where to emit the neutrons from
    const G4double sourceH = 2 * inch;
    const G4double sourceR = 0.2 * inch;
    G4Tubs* const solidSource
        = new G4Tubs("NeutronSource",                   // the name
                     0.0,                               // inner radius
                     sourceR,                           // outer radius
                     sourceH * 0.5,                     // half height
                     0.0 * deg,                         // start angle
                     360.0 * deg);                      // end angle
    G4LogicalVolume* const logicalSource
        = new G4LogicalVolume(solidSource,              // the solid volume
                              fStainlessSteel,          // the material
                              solidSource->GetName());  // the name
    // Place in the exact center of the water tank
    position.set(0.0, 0.0, 0.0);
    new G4PVPlacement(NULL,                             // no rotation
                      position,                         // shift to origin
                      logicalSource,                    // the logical volume
                      logicalSource->GetName(),         // the name
                      logicalTankH2O,                   // the mother volume
                      false,                            // no boolean ops
                      0,                                // copy number
                      overlapChecking);                 // check for overlaps
    
    //
    // Detector Tower
    //
    const G4double polyS = 3.0 * inch;
    const G4double polyH = 18.0 * inch;
    G4Box* const solidPoly
        = new G4Box("Poly",                             // the name
                    polyS,                              // x size
                    polyS,                              // y size
                    polyH);                             // z size
    G4LogicalVolume* const logicalPoly
        = new G4LogicalVolume(solidPoly,                // the solid volume
                              fPolyethylene,            // the material
                              solidPoly->GetName());    // the name
    // The polyethylene detector tower goes just outside the tank at 45 deg
    G4double radiusToPolyCenter = (tankOR / std::sqrt(2.0)) + std::sqrt(2.0) * polyS;
    position.set(-radiusToPolyCenter, radiusToPolyCenter, polyH);
    new G4PVPlacement(NULL,                             // no rotation
                      position,                         // position
                      logicalPoly,                      // the logical volume
                      logicalPoly->GetName(),           // the name
                      logicalWorld,                     // the mother volume
                      false,                            // no boolean ops
                      0,                                // copy number
                      overlapChecking);                 // check for overlaps
    // Create the detector shell
    G4double shellR = 0.3 * inch;
    G4double shellH = 6.5 * inch;
    G4Tubs* const solidShell
        = new G4Tubs("Detector_Shell",                  // the name
                     0.0,                               // inner radius
                     shellR,                            // outer radius
                     shellH * 0.5,                      // half height
                     0.0 * deg,                         // start angle
                     360.0 * deg);                      // end angle
    G4LogicalVolume* const logicalShell
        = new G4LogicalVolume(solidShell,               // the solid volume
                              fStainlessSteel,          // the material
                              solidShell->GetName());   // the name
    // Place in the exact center of the polyethylene tower
    position.set(0.0, 0.0, 0.0);
    new G4PVPlacement(NULL,                             // no rotation
                      position,                         // shift to origin
                      logicalShell,                     // the logical volume
                      logicalShell->GetName(),          // the name
                      logicalPoly,                      // the mother volume
                      false,                            // no boolean ops
                      0,                                // copy number
                      overlapChecking);                 // check for overlaps
    // Create the BF3 detector
    G4double BF3R = 0.2 * inch;
    G4double BF3H = 6.0 * inch;
    G4Tubs* const solidBF3
        = new G4Tubs("Detector_BF3_Core",               // the name
                     0.0,                               // inner radius
                     BF3R,                              // outer radius
                     BF3H * 0.5,                        // half height
                     0.0 * deg,                         // start angle
                     360.0 * deg);                      // end angle
    G4LogicalVolume* const logicalBF3
        = new G4LogicalVolume(solidBF3,                 // the solid volume
                              fBF3_96E,                 // the material
                              solidBF3->GetName());     // the name
    // Place in the exact center of the shell
    position.set(0.0, 0.0, 0.0);
    new G4PVPlacement(NULL,                             // no rotation
                      position,                         // shift to origin
                      logicalBF3,                       // the logical volume
                      logicalBF3->GetName(),            // the name
                      logicalShell,                     // the mother volume
                      false,                            // no boolean ops
                      0,                                // copy number
                      overlapChecking);                 // check for overlaps
    
    return physicalWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FFDetectorConstruction::
DefineMaterials(void)
{
    static G4NistManager* const nist = G4NistManager::Instance();
    
    fAir = nist->FindOrBuildMaterial("G4_AIR");
    fAluminum = nist->FindOrBuildMaterial("G4_Al");
    fGraphite = nist->FindOrBuildMaterial("G4_GRAPHITE");
    fPolyethylene = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
    fStainlessSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    fWater = nist->FindOrBuildMaterial("G4_WATER");
    
    /*// List available materials
    std::vector< G4String > materials = nist->GetNistMaterialNames();
    for(unsigned int i = 0;
        i < materials.size();
        ++i)
    {
        G4cout << materials[i] << G4endl;
    }*/
    
    //
    // Define the 20% enriched UO2
    //
    // First we need to start by creating the isotopes
    G4double const U235Enrichment = 0.2;
    G4double const U238Enrichment = 0.8;
    G4Isotope* const iU235
        = new G4Isotope("iU235",                        // name
                        92,                             // ZZZ
                        235,                            // AAA
                        235.053930 * (g / mole));       // molecular weight
    G4Isotope* const iU238
        = new G4Isotope("iU238",                        // name
                        92,                             // ZZZ
                        238,                            // AAA
                        238.050788 * (g / mole));       // molecular weight
    // Now create the elements and add the isotopes
    G4Element* const U235
        = new G4Element("U235",                         // name
                        "U235",                         // symbol
                        1);                             // number of isotopes
    U235->AddIsotope(iU235,                             // isotope
                     1.0);                              // abundance
    G4Element* const U238
        = new G4Element("U238",                         // name
                        "U238",                         // symbol
                        1);                             // number of isotopes
    U238->AddIsotope(iU238,                             // isotope
                     1.0);                              // abundance
    G4Element* const oxygen = nist->FindOrBuildElement("O");
    // Calculate the mass fractions
    const G4double UO2MolecularWeight = U235->GetA() * U235Enrichment
                                        + U238->GetA() * U238Enrichment
                                        + oxygen->GetA() * 2;
    const G4double U235MassFraction = (U235->GetA() * U235Enrichment)
                                      / UO2MolecularWeight;
    const G4double U238MassFraction = (U238->GetA() * U238Enrichment)
                                      / UO2MolecularWeight;
    const G4double oxygenMassFraction = (oxygen->GetA() * 2)
                                        / UO2MolecularWeight;
    // create the material and add the elements
    fUO2_20E = new G4Material("UO2_20E",                // name
                              10.97 * (g / cm3),        // density
                              3);                       // number of components
    fUO2_20E->AddElement(U235,                          // element
                         U235MassFraction);             // mass fraction
    fUO2_20E->AddElement(U238,                          // element
                         U238MassFraction);             // mass fraction
    fUO2_20E->AddElement(oxygen,                        // element
                         oxygenMassFraction);           // mass fraction
    
    //
    // Define the BF3
    //
    // The BF3 is B-10 enriched
    // http://www.orau.org/ptp/collection/proportional%20counters/bf3info.htm
    G4double const B10Enrichment = 0.96;
    G4double const B11Enrichment = 0.04;
    G4Isotope* const iB10
        = new G4Isotope("iB10",                         // name
                        5,                              // ZZZ
                        10,                             // AAA
                        10.0129370 * (g / mole));       // molecular weight
    G4Isotope* const iB11
        = new G4Isotope("iB11",                         // name
                        5,                              // ZZZ
                        11,                             // AAA
                        11.0093054 * (g / mole));       // molecular weight
    // Now create the elements and add the isotopes
    G4Element* const B10
        = new G4Element("B10",                          // name
                        "B10",                          // symbol
                        1);                             // number of isotopes
    B10->AddIsotope(iB10,                               // isotope
                     1.0);                              // abundance
    G4Element* const B11
        = new G4Element("B11",                          // name
                        "B11",                          // symbol
                        1);                             // number of isotopes
    B11->AddIsotope(iB11,                               // isotope
                     1.0);                              // abundance
    G4Element* const flouride = nist->FindOrBuildElement("F");
    // Calculate the mass fractions
    const G4double BF3MolecularWeight = B10->GetA() * B10Enrichment
                                        + B11->GetA() * B11Enrichment
                                        + flouride->GetA() * 3;
    const G4double B10MassFraction = (B10->GetA() * B10Enrichment)
                                     / BF3MolecularWeight;
    const G4double B11MassFraction = (B11->GetA() * B11Enrichment)
                                     / BF3MolecularWeight;
    const G4double flourideMassFraction = (flouride->GetA() * 3)
                                          / BF3MolecularWeight;
    // create the material and add the elements
    fBF3_96E = new G4Material("BF3_96E",                // name
                              2.5 * (kg / m3),          // density
                              3);                       // number of components
    fBF3_96E->AddElement(B10,                           // element
                         B10MassFraction);              // mass fraction
    fBF3_96E->AddElement(B11,                           // element
                         B11MassFraction);              // mass fraction
    fBF3_96E->AddElement(flouride,                      // element
                         flourideMassFraction);         // mass fraction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FFDetectorConstruction::
PlaceFuelPlate(double x,
               double y,
               G4LogicalVolume* const myLogicalVolume,
               G4LogicalVolume* const parentLogicalVolume)
{
    G4ThreeVector position(x, y);
    std::ostringstream copyName;
    copyName << "Plate@Location    X:" << std::setprecision(2) << x / inch << "    Y:" << y / inch; 
    
    new G4PVPlacement(NULL,                 // no rotation
                      position,             // position
                      myLogicalVolume,      // the logical volume
                      copyName.str(),       // the name
                      parentLogicalVolume,  // the mother volume
                      false,                // no boolean ops
                      fCopyNumber++,        // copy number
                      true);                // check for overlaps
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
FFDetectorConstruction::
~FFDetectorConstruction()
{
    // Nothing here
}


