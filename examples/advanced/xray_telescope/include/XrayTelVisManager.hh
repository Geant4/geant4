//  XrayTelVisManager.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayTelVisManager_h
#define XrayTelVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelVisManager: public G4VisManager {

public:

  XrayTelVisManager ();

private:

  void RegisterGraphicsSystems ();

};

#endif
#endif
