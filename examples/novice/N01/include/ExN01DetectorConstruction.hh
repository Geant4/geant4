// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN01DetectorConstruction.hh,v 1.1.10.1 1999/12/07 20:47:21 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#ifndef ExN01DetectorConstruction_H
#define ExN01DetectorConstruction_H 1

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class ExN01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExN01DetectorConstruction();
    ~ExN01DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

};

#endif

