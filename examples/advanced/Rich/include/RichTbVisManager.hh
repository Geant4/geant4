// Rich advanced example for Geant4
// RichTbVisManager.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////

#ifndef RichTbVisManager_h
#define RichTbVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class RichTbVisManager: public G4VisManager {

public:

  RichTbVisManager ();
  virtual ~RichTbVisManager();
private:

  void RegisterGraphicsSystems ();

};

#endif

#endif
