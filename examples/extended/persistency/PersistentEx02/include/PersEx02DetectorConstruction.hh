// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02DetectorConstruction.hh,v 1.3 1999/11/29 18:15:28 morita Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#ifndef PersEx02DetectorConstruction_h
#define PersEx02DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class PersEx02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    PersEx02DetectorConstruction();
    ~PersEx02DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

};

#endif

