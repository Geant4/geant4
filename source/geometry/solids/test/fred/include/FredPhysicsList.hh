//
// FredPhysicsList.hh
//
// Declaration of Physics for fred
//

#ifndef FredPhysicslist_h
#define FredPhysicslist_h

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class FredPhysicsList : public G4VUserPhysicsList
{
	public:
	FredPhysicsList();
	~FredPhysicsList();
	
	protected:
	void ConstructParticle();
	void ConstructProcess();
	void SetCuts();
};

#endif
