// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst01DetectorConstruction.hh,v 1.1 2001-02-08 08:41:42 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst01DetectorConstruction_h
#define Tst01DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst01DetectorConstruction();
    ~Tst01DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
};

#endif

