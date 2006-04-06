// -------------------------------------------------------------------
// $Id: MicrobeamSteppingAction.hh,v 1.1 2006-04-06 15:32:44 sincerti Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamSteppingAction_h
#define MicrobeamSteppingAction_h 1

#include "G4UserSteppingAction.hh"

#include "MicrobeamRunAction.hh"
#include "MicrobeamDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MicrobeamSteppingAction : public G4UserSteppingAction
{
public:
  MicrobeamSteppingAction(MicrobeamRunAction* ,MicrobeamDetectorConstruction*);
  ~MicrobeamSteppingAction();
  
  void UserSteppingAction(const G4Step*);
  
private:
  MicrobeamRunAction*            Run;
  MicrobeamDetectorConstruction* Detector;
  G4float massPhantom;

};

#endif




