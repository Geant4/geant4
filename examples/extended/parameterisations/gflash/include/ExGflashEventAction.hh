#ifndef ExGflashEventAction_h
#define ExGflashEventAction_h

#include "globals.hh"
#include "G4UserEventAction.hh"


class ExGflashEventAction: public G4UserEventAction {
	public:
	ExGflashEventAction();
	~ExGflashEventAction();
	
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);
	private:
	G4int nevent;
	G4double dtime;
	G4int calorimeterCollectionId;
};
#endif
