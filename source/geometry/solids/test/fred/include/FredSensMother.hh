//
// FredSensMother.hh
//
// Definition of Fred's sensitive mother volume
//

#ifndef FredSensMother_hh
#define FredSensMother_hh

#include "G4VSensitiveDetector.hh"
#include "FredHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class FredSensMother : public G4VSensitiveDetector
{
	public:
	FredSensMother( G4String name );
	~FredSensMother();
	
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
