#ifndef hTestTrackingAction_h
#define hTestTrackingAction_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestTrackingAction
//  
// Description: Implementation file for control on MC truth 
//
// Author:      V.Ivanchenko 13/03/01
//
// Modified: 
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestHisto.hh"
#include "G4UserTrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestTrackingAction : public G4UserTrackingAction
{

  public:
    hTestTrackingAction();
    ~hTestTrackingAction();

    void PreUserTrackingAction(const G4Track*);
    void PostUserTrackingAction(const G4Track*) {};

  private:

    hTestHisto* theHisto;

};

#endif

