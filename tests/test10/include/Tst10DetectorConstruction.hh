// $Id: Tst10DetectorConstruction.hh,v 1.1 1999-01-08 16:35:32 gunter Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserDetectorConstruction
//      for constructing all particles and processes.
//
//	History
//        first version              09 Sept. 1998 by S.Magni
// ------------------------------------------------------------

#ifndef Tst10DetectorConstruction_h
#define Tst10DetectorConstruction_h 1

#include "Tst10DetectorMessenger.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst10DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst10DetectorConstruction();
    ~Tst10DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
		 G4VPhysicalVolume* SelectDetector (G4String val);
		 void SwitchDetector (void);
		 void SetMaterial( void );
	private:
	   Tst10DetectorMessenger* detectorMessenger;
		 G4VSolid* aVolume;
		 G4VPhysicalVolume* PhysicalVolume;
		 G4Material* Water;
		 G4Material* Water1;
		 G4OpticalSurface* aSurface;
};

#endif

