// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em4DetectorConstruction.hh,v 1.1 1999-10-12 11:26:55 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em4DetectorConstruction_h
#define Em4DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em4DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Em4DetectorConstruction();
   ~Em4DetectorConstruction();
     
    G4VPhysicalVolume* Construct();
};

#endif

