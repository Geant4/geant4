// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01EventAction.hh,v 1.1 2001-03-27 16:21:26 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F01EventAction_h
#define F01EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class F01RunAction;
class F01EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F01EventAction : public G4UserEventAction
{
  public:
    F01EventAction(F01RunAction* F01RA);
   ~F01EventAction();

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
    void SetTr();
    void SetRef();
    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int val)  {printModulo = val;};
        
  private:
    G4int    calorimeterCollID;
    F01EventActionMessenger*  eventMessenger;
    F01RunAction* runaction;
    G4int verboselevel;
    G4double nstep,nstepCharged,nstepNeutral;
    G4double Nch,Nne;
    G4double NE,NP;
    G4double Transmitted,Reflected ;
    
    G4String drawFlag;
    G4int    printModulo;             
};

#endif

    
