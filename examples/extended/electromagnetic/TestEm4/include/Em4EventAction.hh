// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em4EventAction.hh,v 1.2 1999-12-15 14:49:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em4EventAction_h
#define Em4EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Em4RunAction;
class Em4EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em4EventAction : public G4UserEventAction
{
  public:
    Em4EventAction(Em4RunAction*);
   ~Em4EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void addEdep(G4double Edep)     {TotalEnergyDeposit += Edep;};      
    G4double GetEnergyDeposit()     {return TotalEnergyDeposit;};    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int val)  {printModulo = val;};
        
  private:
    Em4RunAction* Em4Run;            
    G4double TotalEnergyDeposit;   // Energy deposited in c6f6
    G4String drawFlag;             // control the drawing of event
    G4int                     printModulo;          
    Em4EventActionMessenger*  eventMessenger;
};

#endif

    
