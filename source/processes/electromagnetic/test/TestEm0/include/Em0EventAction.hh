// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0EventAction.hh,v 1.3 1999-12-15 14:51:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em0EventAction_h
#define Em0EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Em0EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em0EventAction : public G4UserEventAction
{
  public:
    Em0EventAction();
   ~Em0EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void addEdep(G4double Edep)    {TotalEnergyDeposit += Edep;};      
    G4double GetEnergyDeposit()    {return TotalEnergyDeposit;};    
    void SetDrawFlag(G4String val) {drawFlag = val;};
    
  private: 
    G4double                  TotalEnergyDeposit;   
    G4String                  drawFlag;             
    Em0EventActionMessenger*  eventMessenger;
};

#endif

    
