// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISAStackingActionMessenger_h
#define LISAStackingActionMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"

class LISAStackingAction;

class LISAStackingActionMessenger: public G4UImessenger {

  public:
    LISAStackingActionMessenger(LISAStackingAction*);
    ~LISAStackingActionMessenger();

    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    LISAStackingAction* stackingAction;   

    G4UIcmdWithABool* SetPriSurvey;
    G4UIcmdWithABool* SetPartSurvey;

};

#endif
