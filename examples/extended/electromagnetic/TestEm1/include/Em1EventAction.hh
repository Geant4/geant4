// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1EventAction.hh,v 1.2 1999-12-15 14:48:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em1EventAction_h
#define Em1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Em1RunAction;
class Em1EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em1EventAction : public G4UserEventAction
{
  public:
    Em1EventAction(Em1RunAction*);
   ~Em1EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void addEdep(G4double Edep)    {TotalEnergyDeposit += Edep;};      
    G4double GetEnergyDeposit()    {return TotalEnergyDeposit;};    
    void SetDrawFlag(G4String val) {drawFlag = val;};
    void SetPrintModulo(G4int val) {printModulo = val;};
            
    
  private:
    Em1RunAction*             Em1Run;
    G4double                  TotalEnergyDeposit;   
    G4String                  drawFlag;
    G4int                     printModulo;                    
    Em1EventActionMessenger*  eventMessenger;
};

#endif

    
