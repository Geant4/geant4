//
// FredVisManager
//
// Definition of visualization manager for Fred
//

#ifndef FredVisManager_H
#define FredVisManager_H

#include "G4VisManager.hh"

class FredVisManager : public G4VisManager
{
	public:
	FredVisManager() {}
	~FredVisManager() {}
	
	private:
	void RegisterGraphicsSystems();
}; 

#endif
