// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F03SteppingAction.hh,v 1.1 2001-06-08 11:55:40 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F03SteppingAction_h
#define F03SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class F03DetectorConstruction;
class F03RunAction;
class F03EventAction;
class F03SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F03SteppingAction : public G4UserSteppingAction
{
  public:
    F03SteppingAction(F03DetectorConstruction*, F03EventAction*,
                      F03RunAction* );
   ~F03SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    F03DetectorConstruction* detector;
    F03EventAction*          eventaction;
    F03RunAction*            runaction;
    F03SteppingMessenger*    steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif
