// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISASteppingActionMessenger class                                *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************

#include "LISASteppingActionMessenger.hh"

#include "LISASteppingAction.hh"


LISASteppingActionMessenger::LISASteppingActionMessenger
  (LISASteppingAction* stepAct) : steppingAction(stepAct){


  newDir = new G4UIdirectory("/surveys/");
  newDir->SetGuidance("Particle survey control commands.");

  // set flag for spectrum saving
  SetFlagSpectrum = new G4UIcmdWithABool("/surveys/surveyTestMasses",this);
  SetFlagSpectrum->SetGuidance("Dump energies (eV) for e- and hadrons");
  SetFlagSpectrum->SetGuidance("entering and leaving the test masses.");
  SetFlagSpectrum->SetGuidance("'electrons_in.out' : e- entering TM");
  SetFlagSpectrum->SetGuidance("'electrons_out.out': e- leaving TM");
  SetFlagSpectrum->SetGuidance("'hadrons_in.out'   : hadrons entering TM");
  SetFlagSpectrum->SetGuidance("'hadrons_out.out'  : hadrons leaving TM");
  SetFlagSpectrum->SetGuidance("Default = false");
  SetFlagSpectrum->SetParameterName("FlagSpectrum", false);
  SetFlagSpectrum->AvailableForStates(G4State_Idle);

}


LISASteppingActionMessenger::~LISASteppingActionMessenger() {

  delete SetFlagSpectrum;  
  delete newDir;

}


void LISASteppingActionMessenger::SetNewValue
   (G4UIcommand* command, G4String newValue) { 

  if(command == SetFlagSpectrum) {
    G4int vl;
    const char* t = newValue;
    std::istrstream is((char*)t);
    is >> vl;
    steppingAction->SetFlagSpectrum(vl!=0);
  }



}

