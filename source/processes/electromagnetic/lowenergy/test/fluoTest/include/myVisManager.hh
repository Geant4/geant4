//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Example Visualization Manager implementing virtual function
//   RegisterGraphicsSystems.  Exploits C-pre-processor variables
//   G4VIS_USE_DAWN, etc., which are set by the GNUmakefiles if
//   environment variables of the same name are set.

// So all you have to do is set environment variables and compile and
//   instantiate this in your main().

// Alternatively, you can implement an empty function here and just
//   register the systems you want in your main(), e.g.:
//   G4VisManager* myVisManager = new MyVisManager;
//   myVisManager -> RegisterGraphicsSystem (new MyGraphicsSystem);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef myVisManager_h
#define myVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myVisManager: public G4VisManager {

public:

  myVisManager ();

private:

  void RegisterGraphicsSystems ();

};

#endif

#endif
