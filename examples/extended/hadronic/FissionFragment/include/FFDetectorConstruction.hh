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
    G4Material* air;
    G4Material* aluminum;
    G4Material* BF3_96E;
    G4Material* concrete;
    G4Material* stainlessSteel;
    G4Material* polyethylene;
    G4Material* UO2_20E;
    G4Material* water;

// Functions
    void DefineMaterials(void);
};

#endif //FFDETECTORCONSTRUCTION


