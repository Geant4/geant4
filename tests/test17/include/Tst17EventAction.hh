// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17EventAction.hh,v 1.2 1999-12-15 14:54:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst17EventAction_h
#define Tst17EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Tst17RunAction;
class Tst17EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst17EventAction : public G4UserEventAction
{
  public:
    Tst17EventAction(Tst17RunAction* Tst17RA);
   ~Tst17EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
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

    G4int    calorimeterCollID;                // Hits collection ID
    Tst17EventActionMessenger*  eventMessenger;
    Tst17RunAction* runaction;            // pointer to the run action
    G4int verboselevel;
    G4double nstep,nstepCharged,nstepNeutral ;      // #steps in the event
    G4double Nch,Nne ;
    G4double NE,NP;
    G4String drawFlag;
    G4double Transmitted,Reflected ;
};

#endif

    







