//
// FredEventAction.hh
//
// Definition of Fred's user action
//

#ifndef FredEventAction_H
#define FredEventAction_H

#include "G4UserEventAction.hh"
#include "globals.hh"

class FredMessenger;

class FredEventAction : public G4UserEventAction
{
	public:
	FredEventAction( FredMessenger *messenger );
	~FredEventAction();
	
	public:
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);
	
	private:
	void DrawNormal();
	void DrawShadow();
	
	private:
	G4int hitID, hitIDMother;
	FredMessenger *messenger;
};

#endif
