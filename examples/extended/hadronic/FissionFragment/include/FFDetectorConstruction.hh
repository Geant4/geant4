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
//! \file       FFDetectorConstruction.hh
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Definition of the FFDetectorConstruction class
//!
//  ================ End Documentation Comments ================
//
//  Modified:
//
//  23-06-14                                              BWendt
//  Added method "PlaceFuelPlates"
//
// -------------------------------------------------------------

#ifndef FFDETECTORCONSTRUCTION
#define FFDETECTORCONSTRUCTION

#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"


class FFDetectorConstruction
:   public G4VUserDetectorConstruction
{
public:
// Constructor
    FFDetectorConstruction();
    
// Functions
    virtual G4VPhysicalVolume* Construct();
    
// Destructor
    virtual ~FFDetectorConstruction();

private:
// Fields
    G4Material* fAir;
    G4Material* fAluminum;
    G4Material* fBF3_96E;
    G4Material* fGraphite;
    G4Material* fStainlessSteel;
    G4Material* fPolyethylene;
    G4Material* fUO2_20E;
    G4Material* fWater;
    unsigned int fCopyNumber;

// Functions
    void DefineMaterials(void);
    void PlaceFuelPlate(double x,
                        double y,
                        G4LogicalVolume* const myLogicalVolume,
                        G4LogicalVolume* const parentLogicalVolume);
};

#endif //FFDETECTORCONSTRUCTION


