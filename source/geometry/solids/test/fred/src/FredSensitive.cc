//
// FredSensitive.cc
//
// Implementation of our sensitive test volume
//

#include "FredSensitive.hh"
#include "G4StepStatus.hh"

//
// Constructor
//
FredSensitive::FredSensitive( G4String name ) : G4VSensitiveDetector( name )
{
	//
	// Now we play the c++ guessing game: which inherited
	// parameters do we fill? Should be obvious, but,
	// the only clue I got was in the examples.
	//
	// IMHO: it would be
	// better to have an inherited class with more arguments
	// whose constructor is explicitly invoked here.
	// Then there would be no ambiguities.
	//
	collectionName.insert( "fredsStuff" );
}

//
// Destructor
//
// Odd we don't bother deleting the hits here... shouldn't we?
// Perhaps it doesn't matter, since G4 is pretty much static when
// it comes to detector configuration. But, we're all c++ experts,
// so we need to go through the moves.
//
FredSensitive::~FredSensitive() {;}

//
// Initialize
//
// Despite the misnomer, this is done at the beginning of each event.
// HCE is a meta-list of hit collections where one is allowed to store
// stuff.
//
void FredSensitive::Initialize( G4HCofThisEvent *HCE )
{
	hits = new FredHitCollection( SensitiveDetectorName, collectionName[0] );

	static int HCID = -1;
	if (HCID < 0) HCID = GetCollectionID(0);
	HCE->AddHitsCollection( HCID, hits );
}

//
// ProcessHits
//
// Another misnomer: this member function processes one step in the 
// sensitive volume.
//
G4bool FredSensitive::ProcessHits( G4Step *step, G4TouchableHistory *notUsed )
{
	G4bool imHit = false;
	
	//
	// Okay, for test purposes, I want to store hits that intersect
	// the boundary of the detector
	//
	if (step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
		FredHit *hit = new FredHit( step->GetPreStepPoint()->GetPosition(), true, step->GetTrack() );
		hits->insert( hit );
		imHit = true;
	}
	if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
		FredHit *hit = new FredHit( step->GetPostStepPoint()->GetPosition(), false, step->GetTrack() );
		hits->insert( hit );
		imHit = true;
	}
	
	return imHit;
}

