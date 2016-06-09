//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
