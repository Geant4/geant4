//
// FredTrackCheck.hh
//
// Implementation of Fred's track hit checker
//

#ifndef FredTrackCheck_hh
#define FredTrackCheck_hh

#include "globals.hh"
#include "G4Track.hh"

#include "g4rw/tpsrtvec.h"


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
	G4RWTPtrSortedVector<FredTrackData> hitList;
};

#endif
