// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14EventAction.hh,v 1.1 1999-05-29 14:12:04 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst14EventAction_h
#define Tst14EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Tst14RunAction;
class Tst14EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst14EventAction : public G4UserEventAction
{
  public:
    Tst14EventAction(Tst14RunAction* Tst14RA);
   ~Tst14EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    G4int GetEventno();
    void setEventVerbose(G4int level);
    
    void CountStepsCharged() ;
    void CountStepsNeutral() ;
    void AddCharged() ;
    void AddNeutral() ;
    void AddE();
    void AddP();   
    void SetTr() ;
    void SetRef() ;
  private:

    G4int calorimeterCollID;                // Hits collection ID
    Tst14EventActionMessenger*  eventMessenger;
    Tst14RunAction* runaction;            // pointer to the run action
    G4int verboselevel;
    G4double nstep,nstepCharged,nstepNeutral ;      // #steps in the event
    G4double Nch,Nne ;
    G4double NE,NP;
    G4double Transmitted,Reflected ;
};

#endif

    
