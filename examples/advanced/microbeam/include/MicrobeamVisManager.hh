// -------------------------------------------------------------------
// $Id: MicrobeamVisManager.hh,v 1.1 2006-04-06 15:32:44 sincerti Exp $
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
