//
//
//    ****************************************
//    *                                      *
//    *    ThyroidDetectorConstruction.hh    *
//    *                                      *
//    ****************************************
//
//

#ifndef ThyroidDetectorConstruction_h
#define ThyroidDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4EllipticalTube.hh"
class G4VPhysicalVolume;
class ThyroidSD;
class G4EllipticalTube;

class ThyroidDetectorConstruction : public G4VUserDetectorConstruction
{
  
 public:
	ThyroidDetectorConstruction(); ~ThyroidDetectorConstruction();
  
 public:
	const G4double m_theDx;
	const G4double m_theDy;
        const G4double m_theDz;

	G4String m_SDName;
        G4int m_pVoxelID;

 public:
	G4VPhysicalVolume* Construct();
};
#endif

