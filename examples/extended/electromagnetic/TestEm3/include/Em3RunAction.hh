// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3RunAction.hh,v 1.4 2000-12-07 12:38:24 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em3RunAction_h
#define Em3RunAction_h 1

#include "Em3DetectorConstruction.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;
class Em3RunActionMessenger;

class Em3RunAction : public G4UserRunAction
{
  public:
  
    Em3RunAction(Em3DetectorConstruction*);
   ~Em3RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);
    
    void fillPerEvent(G4int,G4double,G4double);
    
    void  SetRndmFreq(G4int   val)  {saveRndm = val;}
    G4int GetRndmFreq()             {return saveRndm;}
                
  private:
    void bookHisto();
    
  private:
    
    G4double sumEAbs[MaxAbsor], sum2EAbs[MaxAbsor];
    G4double sumLAbs[MaxAbsor], sum2LAbs[MaxAbsor];

    Em3DetectorConstruction* Detector;    
    Em3RunActionMessenger*   runMessenger;        
    G4int saveRndm;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
void Em3RunAction::fillPerEvent(G4int kAbs, G4double EAbs, G4double LAbs)
{      
  //accumulate statistic
  //
  sumEAbs[kAbs] += EAbs; sum2EAbs[kAbs] += EAbs*EAbs;
  sumLAbs[kAbs] += LAbs; sum2LAbs[kAbs] += LAbs*LAbs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

