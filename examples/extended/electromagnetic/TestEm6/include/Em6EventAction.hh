// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// Event action definition.
// Class Description - end
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em6EventAction_h
#define Em6EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Em6RunAction;
class Em6EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em6EventAction : public G4UserEventAction
{
public: // Without description

    Em6EventAction(Em6RunAction* Em6RA);
   ~Em6EventAction();

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
    void SetTr();
    void SetRef();
    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    
  private:
    G4int    calorimeterCollID;
    Em6EventActionMessenger*  eventMessenger;
    Em6RunAction* runaction;
    G4int verboselevel;
    G4double nstep,nstepCharged,nstepNeutral;
    G4double Nch,Nne;
    G4double NE,NP;
    G4double Transmitted,Reflected ;
    
    G4String drawFlag;
};

#endif

    
