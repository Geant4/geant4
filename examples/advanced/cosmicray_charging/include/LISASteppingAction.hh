// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISASteppingAction_h
#define LISASteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class LISAEventAction;
class LISASteppingActionMessenger;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class LISASteppingAction : public G4UserSteppingAction {

  public:
    LISASteppingAction();
    virtual ~LISASteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  public:
    inline void Discharge(void)          {charge_in[0]=0; charge_out[0]=0;
                                         charge_in[1]=0; charge_out[1]=0;};

    inline G4int GetChargeIn(const G4int tm)       {return charge_in[tm];}
    inline G4int GetChargeOut(const G4int tm)      {return charge_out[tm];}


    inline void SetFlagSpectrum(const G4bool val)  {FlagSpectrum  = val;}

  private:
    G4int charge_in[2];
    G4int charge_out[2];

    LISASteppingActionMessenger* steppingMessenger;
    G4bool FlagSpectrum;



};

#endif
