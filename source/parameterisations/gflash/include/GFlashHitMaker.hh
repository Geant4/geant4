// Created by Joanna Weng 9.11.2004 

#ifndef GFlashHitMaker_h
#define GFlashHitMaker_h 1
// has singleton semantics

#include "G4TouchableHandle.hh"
#include "G4Navigator.hh"

#include "GFlashEnergySpot.hh"
#include "G4GFlashSpot.hh"
#include "G4FastTrack.hh"

class GFlashHitMaker 
{
	public:
	GFlashHitMaker();
	~GFlashHitMaker();
	
	void make(GFlashEnergySpot * aSpot, const G4FastTrack * aT );
	
	private:  
	G4TouchableHandle fTouchableHandle;
	G4Navigator *fpNavigator;
	G4bool fNaviSetup;
	
	private:
	GFlashHitMaker(const GFlashHitMaker & )
	{
		G4Exception("GFlashHitMaker copy constructor not publicly available");
	}
	GFlashHitMaker & operator = (const GFlashHitMaker & )
	{
		G4Exception("GFlashHitMaker asignment operator not publicly available");
		return *this;
	}
	
};
#endif

