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

#ifndef Test17EventAction_h
#define Test17EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Test17RunAction;
class Test17EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17EventAction : public G4UserEventAction
{
public: // Without description

    Test17EventAction(Test17RunAction* Test17RA);
   ~Test17EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    G4int GetEventno();
    void setEventVerbose(G4int level);
    
    void CountStepsCharged() ;
    void CountStepsNeutral() ;
    void AddCharged() ;
    void AddNeutral() ;
    void AddE(G4double En);
    void AddP();   
    void SetTr();
    void SetRef();
    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    
  private:
    G4int    calorimeterCollID;
    Test17EventActionMessenger*  eventMessenger;
    Test17RunAction* runaction;
    G4int verboselevel;
    G4double nstep,nstepCharged,nstepNeutral;
    G4double Nch,Nne;
    G4double NE,NP;
    G4double Transmitted,Reflected ;
    G4double EnergyDeposition ;

    G4String drawFlag;
};

#endif

    
