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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISADetectorConstruction_h
#define LISADetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "G4UserLimits.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Region.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "globals.hh"

#include <cmath>


class LISADetectorConstruction : public G4VUserDetectorConstruction {

  public:
    LISADetectorConstruction();
    ~LISADetectorConstruction();

    G4VPhysicalVolume* Construct();

  private:
    void ConstructMaterials();
    G4VPhysicalVolume* ConstructDetector();


  private:
  // pointers to materials
  G4Material *vacuum;
  G4Material *Al6061;
  G4Material *AlHoneycomb;
  G4Material *Scell;
  G4Material *MLImat;
  G4Material *molybdenum;
  G4Material *TiAlloy;
  G4Material *gold;
  G4Material *AuPt;
  G4Material *CFRP;
  G4Material *ULEglass;
  G4Material *SHAPAL;
  G4Material *SiC;
  G4Material *foam;
  

};

#endif

