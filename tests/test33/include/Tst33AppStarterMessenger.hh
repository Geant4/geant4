#ifndef Tst33AppStarterMessenger_hh
#define Tst33AppStarterMessenger_hh Tst33AppStarterMessenger_hh

#include "g4std/map"
#include "G4UImessenger.hh"

class G4UIcmdWithAnInteger;
class Tst33AppStarter;
class Tst33VApplication;

typedef G4std::map<G4UIcmdWithAnInteger *, G4String> Tst33SimComands;

class Tst33AppStarterMessenger : public G4UImessenger{
public:
  explicit Tst33AppStarterMessenger(Tst33AppStarter &);
  virtual ~Tst33AppStarterMessenger();
  virtual void SetNewValue(G4UIcommand* pCmd,G4String szValue);

private:
  Tst33AppStarterMessenger(const Tst33AppStarterMessenger &);
  Tst33AppStarterMessenger &operator=(const Tst33AppStarterMessenger &);
  Tst33AppStarter &fAppStarter;
  G4UIcommand *fMassGeoCmd;
  G4UIcommand *fParallelGeoCmd;
  G4UIcommand *fScoringCmd;
  G4UIcommand *fImpCmd;
  G4UIcommand *fWWRCmd;
  G4UIcommand *fClearSmaplingCmd;
  G4UIcommand *fConfigureSamplingCmd;
  G4UIcommand *fVisAppComand;
  G4UIcmdWithAnInteger *fTimedAppComand;
  G4UIcommand *fPostRunCmd;
  G4UIcmdWithAnInteger *fRunCmd;
};

#endif
