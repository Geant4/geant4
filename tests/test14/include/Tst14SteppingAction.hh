// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14SteppingAction.hh,v 1.1 1999-05-29 14:12:07 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst14SteppingAction_h
#define Tst14SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class Tst14DetectorConstruction;
class Tst14RunAction;
class Tst14EventAction;
class Tst14SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst14SteppingAction : public G4UserSteppingAction
{ 
  public:
    Tst14SteppingAction(Tst14DetectorConstruction* DET,Tst14EventAction* EA,
                        Tst14RunAction* RA);
   ~Tst14SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    Tst14DetectorConstruction* detector;
    Tst14EventAction* eventaction;
    Tst14RunAction* runaction;
    Tst14SteppingMessenger* steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif

