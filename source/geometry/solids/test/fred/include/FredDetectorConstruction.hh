//
// FredDetectorConstruction.hh
//
// Definition of fred's detector
//

#ifndef FredDetectorConstruction_H
#define FredDetectorConstruction_H

#include "FredMessenger.hh"

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class FredDetectorConstruction : public G4VUserDetectorConstruction
{
	public:
	FredDetectorConstruction( FredMessenger *ourMessenger );
	~FredDetectorConstruction();
	
	G4VPhysicalVolume *Construct();
	
	inline G4VSolid	  *GetTestVolume() {return testVolume;}
	
	private:
	FredMessenger	*messenger;
	G4VSolid	*testVolume;
};

#endif
