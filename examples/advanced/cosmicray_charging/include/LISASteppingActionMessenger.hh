// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISASteppingActionMessenger_h
#define LISASteppingActionMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"
#include "globals.hh"

class LISASteppingAction;

class LISASteppingActionMessenger: public G4UImessenger {

  public:
    LISASteppingActionMessenger(LISASteppingAction*);
    ~LISASteppingActionMessenger();

    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    LISASteppingAction* steppingAction;   

    G4UIcmdWithABool* SetFlagSpectrum;

    G4UIdirectory* newDir;

};

#endif
