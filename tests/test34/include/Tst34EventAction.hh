#ifndef Tst34EventAction_h
#define Tst34EventAction_h

#include "globals.hh"
#include "G4UserEventAction.hh"


class Tst34EventAction: public G4UserEventAction {
	public:
	Tst34EventAction();
	~Tst34EventAction();
	
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);
	private:
	G4int nevent;
	G4double dtime;
	G4int calorimeterCollectionId;
};
#endif
