#ifndef hTestEventAction_h
#define hTestEventAction_h 1
 
// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestEventAction -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4UserEventAction.hh"
#include "hTestDetectorConstruction.hh"
#include "G4Event.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestEventAction : public G4UserEventAction
{
public: // Without description

    hTestEventAction(const hTestDetectorConstruction*);
   ~hTestEventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    inline void SetDrawFlag(G4String val)  {drawFlag = val;};
    inline G4int GetEventNo() const {return nEvt;};
    
  private:

    const hTestDetectorConstruction* theDet;
    G4int nEvt;
    G4int verbose;
    G4String drawFlag;
};

#endif

    
