//
// FredSensitive
//
// Definition of Fred's test sensitive detector
//

#ifndef FredSensitive_hh
#define FredSensitive_hh

#include "G4VSensitiveDetector.hh"
#include "FredHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class FredSensitive : public G4VSensitiveDetector
{
	public:
	FredSensitive( G4String name );
	~FredSensitive();
	
	//
	// G4V virtual functions all of which must be overridden
	//
	void Initialize( G4HCofThisEvent *HCE );
	G4bool ProcessHits( G4Step *step, G4TouchableHistory *hist );

	//
	// We hold onto the hits here
	//
	private:
	FredHitCollection *hits;
};

#endif
