// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em10SteppingAction.hh,v 1.1 2000-07-14 15:51:17 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em10SteppingAction_h
#define Em10SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class Em10DetectorConstruction;
class Em10RunAction;
class Em10EventAction;
class Em10SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em10SteppingAction : public G4UserSteppingAction
{
  public:
    Em10SteppingAction(Em10DetectorConstruction*, Em10EventAction*,
                      Em10RunAction* );
   ~Em10SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    Em10DetectorConstruction* detector;
    Em10EventAction*          eventaction;
    Em10RunAction*            runaction;
    Em10SteppingMessenger*    steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif
