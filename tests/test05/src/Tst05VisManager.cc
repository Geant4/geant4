// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst05VisManager.cc,v 1.2 2000-05-12 14:45:40 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4VIS_USE

#include "Tst05VisManager.hh"

// Supported drivers...

#ifdef G4VIS_USE_RAYTRACER
#include "G4RayTracer.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst05VisManager::Tst05VisManager () {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst05VisManager::RegisterGraphicsSystems () {

#ifdef G4VIS_USE_RAYTRACER
  RegisterGraphicsSystem (new G4RayTracer);
#endif

  if (fVerbose > 0) {
    G4cout <<
      "\nYou have successfully chosen to use the following graphics systems."
	 << G4endl;
    PrintAvailableGraphicsSystems ();
  }
}

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
