// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISAEventAction_h
#define LISAEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <strstream>

class LISAPrimaryGeneratorAction;
class LISASteppingAction;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class LISAEventAction : public G4UserEventAction {

  public:
    LISAEventAction(LISAPrimaryGeneratorAction*, LISASteppingAction*);
    virtual ~LISAEventAction();
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    LISAPrimaryGeneratorAction* genAction;
    LISASteppingAction* stepAction;

  public:
  inline void SetFilename(const G4String& name)  {filename = name;};


  private:
    G4int event_id;
    G4double energy_pri;
    G4int  charge_in[2];
    G4int charge_out[2];
    G4int charge_tot[2];
    const long* seeds;
    G4String filename;

};

#endif


