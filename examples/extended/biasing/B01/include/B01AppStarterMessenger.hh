#ifndef B01AppStarterMessenger_hh
#define B01AppStarterMessenger_hh B01AppStarterMessenger_hh

#include "g4std/map"
#include "B01SimulationFactory.hh"
#include "G4UImessenger.hh"

class G4UIcmdWithAnInteger;
class B01AppStarter;
class B01VSimulation;
class B01VApplication;

typedef G4std::map<G4UIcmdWithAnInteger *, G4String> B01SimComands;

class B01AppStarterMessenger : public G4UImessenger{
public:
  explicit B01AppStarterMessenger(B01AppStarter *);
  virtual ~B01AppStarterMessenger();
  virtual void SetNewValue(G4UIcommand* pCmd,G4String szValue);

private:
  B01AppStarterMessenger(const B01AppStarterMessenger &);
  B01AppStarterMessenger &operator=(const B01AppStarterMessenger &);
  B01AppStarter *fAppStarter;
  G4bool fWeightRoulette;
  B01VSimulation *fSim;
  G4bool fApp;
  B01SimulationFactory fSimfac;
  B01SimComands fSimComands;
  G4UIcommand *fVisAppComand;
  G4UIcmdWithAnInteger *fTimedAppComand;
  
};

#endif
