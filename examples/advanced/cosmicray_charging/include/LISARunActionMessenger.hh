// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISARunActionMessenger_h
#define LISARunActionMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include <strstream>
#include "globals.hh"

class LISARunAction;

class LISARunActionMessenger: public G4UImessenger {

  public:
    LISARunActionMessenger(LISARunAction*);
   ~LISARunActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
  private:
    LISARunAction* runAction;   

    G4UIcmdWithABool* SetAutoSeedCmd;

};

#endif
