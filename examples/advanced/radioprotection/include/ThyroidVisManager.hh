//
// $Id: ThyroidVisManager.hh,v 1.2 2003-05-23 11:56:00 francy Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 


#ifndef ThyroidVisManager_h
#define ThyroidVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ThyroidVisManager: public G4VisManager {

public:

  ThyroidVisManager ();

private:

  void RegisterGraphicsSystems ();

};

#endif

#endif
