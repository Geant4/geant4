// -------------------------------------------------------------------
// $Id: MicrobeamVisManager.hh,v 1.3 2006-06-01 22:25:19 sincerti Exp $
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
