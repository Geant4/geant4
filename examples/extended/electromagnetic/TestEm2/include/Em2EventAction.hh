// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2EventAction.hh,v 1.2 1999-12-15 14:48:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em2EventAction_h
#define Em2EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Em2RunAction;
class Em2EventActionMessenger;

class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em2EventAction : public G4UserEventAction
{
  public:
  
    Em2EventAction(Em2RunAction*);
   ~Em2EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {drawFlag    = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};
    
  private:
  
    Em2RunAction*             Em2Run;
    G4String                  drawFlag;
    G4int                     printModulo;          
    Em2EventActionMessenger*  eventMessenger;
};

#endif

    
