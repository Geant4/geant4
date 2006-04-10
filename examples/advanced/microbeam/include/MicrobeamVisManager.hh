// -------------------------------------------------------------------
// $Id: MicrobeamVisManager.hh,v 1.2 2006-04-10 14:47:31 sincerti Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamVisManager_h
#define MicrobeamVisManager_h 1

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MicrobeamVisManager: public G4VisManager {

public:
  MicrobeamVisManager ();

private:
  void RegisterGraphicsSystems ();

};

#endif
