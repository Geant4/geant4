//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelVisManager_h
#define GammaRayTelVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelVisManager: public G4VisManager {

public:

  GammaRayTelVisManager ();

private:

  void RegisterGraphicsSystems ();

};

#endif

#endif
