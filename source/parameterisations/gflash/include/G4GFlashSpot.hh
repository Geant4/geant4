#ifndef G4GFlashSpot_h
#define G4GFlashSpot_h

#include "G4ThreeVector.hh"
#include "G4FastTrack.hh"
#include "GFlashEnergySpot.hh"
#include "G4TouchableHandle.hh"

class G4GFlashSpot
{
	public:
	G4GFlashSpot(const GFlashEnergySpot * aSpot, const G4FastTrack * aTrack, G4TouchableHandle aH)
	: theSpot(aSpot), theTrack(aTrack), theHandle(aH) {}
	
	~G4GFlashSpot() {}
	
	const GFlashEnergySpot * GetEnergySpot() const {return theSpot;}
	
	const G4FastTrack * GetOriginatorTrack() const {return theTrack;}
	
	G4TouchableHandle GetTouchableHandle() const {return theHandle;}
	
	G4ThreeVector GetPosition() const {return GetOriginatorTrack()->GetPrimaryTrack()->GetPosition();}
		
	private:
	
	const GFlashEnergySpot * theSpot;
	const G4FastTrack * theTrack;
	G4TouchableHandle theHandle;
};

#endif
