//
// FredTrackCheck.cc
//
// Implemenation of track hit checker
//

#include "FredTrackCheck.hh"

//
// Constructor
//
FredTrackCheck::FredTrackCheck() {;}

//
// Destructor
//
FredTrackCheck::~FredTrackCheck()
{ 
	hitList.clearAndDestroy();
}


//
// AddInHit
//
// Associated a new "in" hit with the track
//
void FredTrackCheck::AddInHit( const G4Track *track )
{
	//
	// Do we already know about this track?
	//
	FredTrackData	*trackData;
	FredTrackData	*target = new FredTrackData( track->GetTrackID() );
	
	if (trackData = hitList.find( target )) {
	
		//
		// Yup: increment appropriately
		//
		trackData->IncrementIn();
		delete target;
	}
	else {
		//
		// Nope: make one
		//
		target->IncrementIn();
		hitList.insert( target );
	}
}

//
// AddOutHit
//
// Associated a new "out" hit with the track
//
void FredTrackCheck::AddOutHit( const G4Track *track )
{
	//
	// Do we already know about this track?
	//
	FredTrackData	*trackData;
	FredTrackData	*target = new FredTrackData( track->GetTrackID() );
	
	if (trackData = hitList.find( target )) {
	
		//
		// Yup: increment appropriately
		//
		trackData->IncrementOut();
		delete target;
	}
	else {
		//
		// Nope: make one
		//
		target->IncrementOut();
		hitList.insert( target );
	}
}

//
// GetTrackStat
//
// Return statistics on specified track
//
void FredTrackCheck::GetTrackStat( const G4Track *track, G4int *nIn, G4int *nOut )
{
	GetTrackStatID( track->GetTrackID(), nIn, nOut );
}

//
// GetTrackStatID
//
// Return statistics on specified track
//
void FredTrackCheck::GetTrackStatID( G4int trackID, G4int *nIn, G4int *nOut )
{
	FredTrackData	*trackData;
	FredTrackData	target( trackID );
	
	if (trackData = hitList.find( &target )) {
		*nIn = trackData->GetNIn();
		*nOut = trackData->GetNOut();
	}
	else {
		*nIn = 0;
		*nOut = 0;
	}
}


//
// ClearAll
//
// Clear all track data
//
void FredTrackCheck::Clear()
{
	hitList.clearAndDestroy();
}


