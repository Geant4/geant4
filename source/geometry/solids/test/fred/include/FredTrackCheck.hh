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
//
// FredTrackCheck.hh
//
// Implementation of Fred's track hit checker
//

#ifndef FredTrackCheck_hh
#define FredTrackCheck_hh

#include "globals.hh"
#include "G4Track.hh"

#include "g4std/vector"


class FredTrackData {
	public:
	FredTrackData( G4int thisId ) { id = thisId; nIn = 0; nOut = 0; }
	~FredTrackData() {;}
	
	G4bool operator==(const FredTrackData &right) const { return id == right.id; }
        G4bool operator<(const FredTrackData &right) const { return id < right.id; }
	G4bool operator<=(const FredTrackData &right) const { return id <= right.id; }
	
	void IncrementIn() {nIn++;}
	void IncrementOut() {nOut++;}
	G4int GetNIn() { return nIn; }
	G4int GetNOut() { return nOut; }
	
	private:
	G4int	id, nIn, nOut;
};


class FredTrackCheck {
	
	public:
	FredTrackCheck();
	~FredTrackCheck();
	
	void AddInHit( const G4Track *track );
	void AddOutHit( const G4Track *track );
	
	void GetTrackStat( const G4Track *track, G4int *nIn, G4int *nOut );
	void GetTrackStatID( G4int trackID, G4int *nIn, G4int *nOut );
	void Clear( );

	private:
	G4std::vector<FredTrackData*> hitList;
};

#endif
