#ifndef ExGflashVisManager_h
#define ExGflashVisManager_h
using namespace std;
#include "G4VisManager.hh"

class ExGflashVisManager: public G4VisManager {
public:
	ExGflashVisManager ();
private:
	void RegisterGraphicsSystems ();
};

#endif
