
#ifndef BrachyVisManager_h
#define BrachyVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class BrachyVisManager: public G4VisManager {

public:

  BrachyVisManager ();

private:

  void RegisterGraphicsSystems ();

};

#endif

#endif



