// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef ExN06DetectorConstruction_h
#define ExN06DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class ExN06DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    ExN06DetectorConstruction();
    ~ExN06DetectorConstruction();

  public:

    G4VPhysicalVolume* Construct();

  private:

    G4double expHall_x;
    G4double expHall_y;
    G4double expHall_z;

    G4double tank_x;
    G4double tank_y;
    G4double tank_z;

    G4double bubble_x;
    G4double bubble_y;
    G4double bubble_z;

};

#endif /*ExN06DetectorConstruction_h*/
