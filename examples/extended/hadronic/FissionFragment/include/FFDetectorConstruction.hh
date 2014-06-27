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


