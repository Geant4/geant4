// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em4RunAction.hh,v 1.4 2000-12-07 13:02:11 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em4RunAction_h
#define Em4RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;
class Em4RunActionMessenger;

class Em4RunAction : public G4UserRunAction
{
  public:
    Em4RunAction();
   ~Em4RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void  SetRndmFreq(G4int    val) {saveRndm = val;}
    G4int GetRndmFreq()             {return saveRndm;}
        
  private:
    void bookHisto();
    
  private:  
    Em4RunActionMessenger* runMessenger;
    G4int saveRndm;
};

#endif

